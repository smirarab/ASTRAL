#include <iostream>
#include <cstring>
#include "phylonet_coalescent_Polytree_PTNative.h"
#include "x86intrin.h"

#define NUM_BITS_IN_WORD 64
#define BATCH_SIZE 32

using namespace std;

__attribute__ ((always_inline)) inline unsigned long long F(unsigned long long (* __restrict__ x)[3]){
	unsigned long long a = x[0][0], b = x[0][1], c = x[0][2], d = x[1][0], e = x[1][1], f = x[1][2],
			g = x[2][0], h = x[2][1], i = x[2][2];
	return a * ( (a + e + i - 3)  * e * i + (a + f + h - 3)  * f * h )
		+ b * ((b + d + i - 3)  * d * i + (b + f + g - 3)  * f * g )
		+ c * ((c + d + h - 3)  * d * h + (c + e + g - 3)  * e * g );
}

__attribute__((vector)) __attribute__ ((always_inline)) inline unsigned long long F(unsigned short a, unsigned short b, unsigned short c, unsigned short d, 
		unsigned short e, unsigned short f, unsigned short g, unsigned short h, unsigned short i){
	return a * ( (-3LL + a + e + i)  * e * i + (-3LL + a + f + h)  * f * h )
		+ b * ((-3LL + b + d + i)  * d * i + (-3LL + b + f + g)  * f * g )
		+ c * ((-3LL + c + d + h)  * d * h + (-3LL + c + e + g)  * e * g );
}

__attribute__((vector)) __attribute__ ((always_inline)) inline unsigned long long Fpoly(unsigned long long sx0, unsigned long long sx1, unsigned long long sx2,
		unsigned long long sxy0, unsigned long long sxy1, unsigned long long sxy2,
		unsigned long long q0, unsigned long long q1, unsigned long long q2){
	return ((sx1 - q1) * (sx2 - q2) - sxy0 + q1 * q2) * q0 * (q0 - 1LL)
		+ ((sx2 - q2) * (sx0 - q0) - sxy1 + q2 * q0) * q1 * (q1 - 1LL)
		+ ((sx0 - q0) * (sx1 - q1) - sxy2 + q0 * q1) * q2 * (q2 - 1LL);
}

struct Polytree{
	int n = 0, m = 0, listSize = 0, queueSize = 0, nWord = 0;
	int* __restrict__ queue = nullptr;
	unsigned long long* __restrict__ treeBits = nullptr;
	
	~Polytree(){
		if (queue) delete queue;
		if (treeBits) delete treeBits;
	}
	
	unsigned long long compute(const unsigned long long* __restrict__ b) const{
		unsigned long long (* __restrict__ lst)[3] = new unsigned long long[listSize][3]{};
		unsigned long long (* __restrict__ stk)[3] = new unsigned long long[n + 1][3]{};
	
		unsigned long long weight = 0;
		int stackEnd = 0, listEnd = n, treeCnt = 0;
		unsigned long long treeTotal[3] = {};
		
		long long checksum = 0;
		for (int i = 0; i < n; i++){
			lst[i][0] = (b[0 * nWord + i / NUM_BITS_IN_WORD] & (1LL << (i % NUM_BITS_IN_WORD))) ? 1 : 0;
			lst[i][1] = (b[1 * nWord + i / NUM_BITS_IN_WORD] & (1LL << (i % NUM_BITS_IN_WORD))) ? 1 : 0;
			lst[i][2] = (b[2 * nWord + i / NUM_BITS_IN_WORD] & (1LL << (i % NUM_BITS_IN_WORD))) ? 1 : 0;
		}
		
		for (int i = 0; i < queueSize; i++){
			int cmd = queue[i];
			if (cmd == -1) {
				for (int j = 0; j < 3; j++){
					treeTotal[j] = 0;
					for (int k = 0; k < nWord; k++){
						treeTotal[j] += __builtin_popcountll(treeBits[treeCnt * nWord + k] & b[j * nWord + k]);
					}
				}
				treeCnt++;
				continue;
			}
			if (cmd & 1){
				int oldStackEnd = stackEnd;
				int numChildren = cmd >> 5;
				unsigned long long tempWeight = 0;
				for (int k = 0; k < 3; k++) stk[oldStackEnd][k] = treeTotal[k];
				for (int j = oldStackEnd - numChildren; j < oldStackEnd; j++){
					for (int k = 0; k < 3; k++) stk[oldStackEnd][k] -= stk[j][k];
				}
				if (cmd & 8){
					if (numChildren == 2){
						tempWeight = F(stk + oldStackEnd - 2);
					}
					else{
						unsigned long long sx0 = 0, sx1 = 0, sx2 = 0, sxy0 = 0, sxy1 = 0, sxy2 = 0;
						for (int j = oldStackEnd - numChildren; j <= oldStackEnd; j++){
							const unsigned long long q0 = stk[j][0], q1 = stk[j][1], q2 = stk[j][2];
							sx0 += q0; sx1 += q1; sx2 += q2;
							sxy0 += q1 * q2; sxy1 += q2 * q0; sxy2 += q0 * q1;
						}
						for (int j = oldStackEnd - numChildren; j <= oldStackEnd; j++){
							const unsigned long long q0 = stk[j][0], q1 = stk[j][1], q2 = stk[j][2];
							tempWeight += ((sx1 - q1) * (sx2 - q2) - sxy0 + q1 * q2) * q0 * (q0 - 1LL)
								+ ((sx2 - q2) * (sx0 - q0) - sxy1 + q2 * q0) * q1 * (q1 - 1LL)
								+ ((sx0 - q0) * (sx1 - q1) - sxy2 + q0 * q1) * q2 * (q2 - 1LL);
						}
					}
					if (cmd & 16) weight += tempWeight * queue[++i];
					else weight += tempWeight;
				}
				stackEnd -= numChildren;
				if (cmd & 2){
					for (int k = 0; k < 3; k++) stk[stackEnd][k] = treeTotal[k] - stk[oldStackEnd][k];
					stackEnd++;
				}
				if (cmd & 4){
					for (int k = 0; k < 3; k++) lst[listEnd][k] = treeTotal[k] - stk[oldStackEnd][k];
					listEnd++;
				}
			}
			else {
				for (int k = 0; k < 3; k++) stk[stackEnd][k] = lst[cmd >> 1][k];
				stackEnd++;
			}
		}
		delete b;
		delete lst;
		delete stk;
		return weight;
	}
	
	void compute(unsigned long long* __restrict__ result, const unsigned long long* __restrict__ b) const{
		__attribute__((aligned(64))) unsigned short (* __restrict__ lst)[3 * BATCH_SIZE] = new __attribute__((aligned(64))) unsigned short[listSize][3 * BATCH_SIZE]{};
		__attribute__((aligned(64))) unsigned short (* __restrict__ stk)[3 * BATCH_SIZE] = new __attribute__((aligned(64))) unsigned short[n + 1][3 * BATCH_SIZE]{};
		
		__attribute__((aligned(64))) unsigned long long weight[BATCH_SIZE] = {};
		int stackEnd = 0, listEnd = n, treeCnt = 0;
		__attribute__((aligned(64))) unsigned short treeTotal[3 * BATCH_SIZE] = {};
		
		for (int i = 0; i < n; i++){
			for (int j = 0; j < 3 * BATCH_SIZE; j++){
				lst[i][j] = (b[j * nWord + i / NUM_BITS_IN_WORD] & (1LL << (i % NUM_BITS_IN_WORD))) ? 1 : 0;
			}
		}
		
		for (int i = 0; i < queueSize; i++){
			int cmd = queue[i];
			if (cmd == -1) {
				for (int j = 0; j < 3 * BATCH_SIZE; j++){
					treeTotal[j] = 0;
					for (int k = 0; k < nWord; k++){
						treeTotal[j] += __builtin_popcountll(treeBits[treeCnt * nWord + k] & b[j * nWord + k]);
					}
				}
				treeCnt++;
				continue;
			}
			if (cmd & 1){
				int oldStackEnd = stackEnd;
				int numChildren = cmd >> 5;
				__attribute__((aligned(64))) unsigned long long tempWeight[BATCH_SIZE] = {};
				#pragma ivdep
				#pragma simd vectorlength(32)
				#pragma vector always
				#pragma vector aligned
				for (int k = 0; k < 3 * BATCH_SIZE; k++) stk[oldStackEnd][k] = treeTotal[k];
				for (int j = oldStackEnd - numChildren; j < oldStackEnd; j++){
					#pragma ivdep
					#pragma simd vectorlength(32)
					#pragma vector always
					#pragma vector aligned
					for (int k = 0; k < 3 * BATCH_SIZE; k++) stk[oldStackEnd][k] -= stk[j][k];
				}
				if (cmd & 8){
					if (numChildren == 2){
						#pragma ivdep
						#pragma simd vectorlength(32)
						#pragma vector always
						#pragma vector aligned
						for (int k = 0; k < BATCH_SIZE; k++){
							tempWeight[k] = F(stk[oldStackEnd - 2][k], stk[oldStackEnd - 2][k + BATCH_SIZE], stk[oldStackEnd - 2][k + 2 * BATCH_SIZE], 
								stk[oldStackEnd - 1][k], stk[oldStackEnd - 1][k + BATCH_SIZE], stk[oldStackEnd - 1][k + 2 * BATCH_SIZE],
								stk[oldStackEnd][k], stk[oldStackEnd][k + BATCH_SIZE], stk[oldStackEnd][k + 2 * BATCH_SIZE]);
						}
					}
					else{
						__attribute__((aligned(64))) unsigned long long sx0[BATCH_SIZE] = {};
						__attribute__((aligned(64))) unsigned long long sx1[BATCH_SIZE] = {};
						__attribute__((aligned(64))) unsigned long long sx2[BATCH_SIZE] = {};
						__attribute__((aligned(64))) unsigned long long sxy0[BATCH_SIZE] = {};
						__attribute__((aligned(64))) unsigned long long sxy1[BATCH_SIZE] = {};
						__attribute__((aligned(64))) unsigned long long sxy2[BATCH_SIZE] = {};
						for (int j = oldStackEnd - numChildren; j <= oldStackEnd; j++){
							#pragma ivdep
							#pragma simd vectorlength(32)
							#pragma vector always
							#pragma vector aligned
							for (int k = 0; k < BATCH_SIZE; k++){
								unsigned long long q0 = stk[j][k], q1 = stk[j][k + BATCH_SIZE], q2 = stk[j][k + 2 * BATCH_SIZE];
								sx0[k] += q0;
								sx1[k] += q1;
								sx2[k] += q2;
								sxy0[k] += q1 * q2;
								sxy1[k] += q2 * q0;
								sxy2[k] += q0 * q1;
							}
						}
						for (int j = oldStackEnd - numChildren; j <= oldStackEnd; j++){
							#pragma ivdep
							#pragma simd vectorlength(32)
							#pragma vector always
							#pragma vector aligned
							for (int k = 0; k < BATCH_SIZE; k++){
								unsigned long long q0 = stk[j][k], q1 = stk[j][k + BATCH_SIZE], q2 = stk[j][k + 2 * BATCH_SIZE];
								tempWeight[k] += Fpoly(sx0[k], sx1[k], sx2[k], sxy0[k], sxy1[k], sxy2[k], q0, q1, q2);
							}
						}
					}
					if (cmd & 16){
						unsigned long long r = queue[++i];
						#pragma ivdep
						#pragma simd vectorlength(32)
						#pragma vector always
						#pragma vector aligned
						for (int k = 0; k < BATCH_SIZE; k++){
							weight[k] += tempWeight[k] * r;
						}
					} 
					else{
						#pragma ivdep
						#pragma simd vectorlength(32)
						#pragma vector always
						#pragma vector aligned
						for (int k = 0; k < BATCH_SIZE; k++){
							weight[k] += tempWeight[k];
						}
					}
				}
				stackEnd -= numChildren;
				if (cmd & 2){
					for (int k = 0; k < 3 * BATCH_SIZE; k++) stk[stackEnd][k] = treeTotal[k] - stk[oldStackEnd][k];
					stackEnd++;
				}
				if (cmd & 4){
					for (int k = 0; k < 3 * BATCH_SIZE; k++) lst[listEnd][k] = treeTotal[k] - stk[oldStackEnd][k];
					listEnd++;
				}
			}
			else {
				for (int k = 0; k < 3 * BATCH_SIZE; k++) stk[stackEnd][k] = lst[cmd >> 1][k];
				stackEnd++;
			}
		}
		delete lst;
		delete stk;
		for (int k = 0; k < BATCH_SIZE; k++){
			result[k] = weight[k];
		}
	}
	
} pt;

JNIEXPORT void JNICALL Java_phylonet_coalescent_Polytree_00024PTNative_cppInit
		(JNIEnv *env, jclass, jint clusterSize, jint listSize, jintArray jQueue, jobjectArray jTreeBits){
	pt.n = clusterSize;
	pt.nWord = (pt.n - 1) / NUM_BITS_IN_WORD + 1;
	pt.listSize = listSize;
	pt.queueSize = env->GetArrayLength(jQueue);
	pt.queue = new int[pt.queueSize]{};
	jint *jptrQueue = env->GetIntArrayElements(jQueue, 0);
	memcpy(pt.queue, jptrQueue, pt.queueSize * sizeof(int));
	env->ReleaseIntArrayElements(jQueue, jptrQueue, 0);
	pt.m = env->GetArrayLength(jTreeBits);
	pt.treeBits = new unsigned long long[pt.m * pt.nWord]{};
	for (int i = 0; i < pt.m; i++){
		jlongArray jBits = (jlongArray) env->GetObjectArrayElement(jTreeBits, i);
		int len = env->GetArrayLength(jBits);
		jlong *jptrBits = env->GetLongArrayElements(jBits, 0);
		memcpy(pt.treeBits + i * pt.nWord, jptrBits, len * sizeof(unsigned long long));
		env->ReleaseLongArrayElements(jBits, jptrBits, 0);
	}
}

JNIEXPORT jlong JNICALL Java_phylonet_coalescent_Polytree_00024PTNative_cppCompute
		(JNIEnv *env, jclass, jlongArray jb1, jlongArray jb2, jlongArray jb3){
	unsigned long long* b = new unsigned long long[3 * pt.nWord]{};
	int len1 = env->GetArrayLength(jb1);
	jlong *jpb1 = env->GetLongArrayElements(jb1, 0);
	memcpy(b + 0 * pt.nWord, jpb1, len1 * sizeof(unsigned long long));
	env->ReleaseLongArrayElements(jb1, jpb1, 0);
	int len2 = env->GetArrayLength(jb2);
	jlong *jpb2 = env->GetLongArrayElements(jb2, 0);
	memcpy(b + 1 * pt.nWord, jpb2, len2 * sizeof(unsigned long long));
	env->ReleaseLongArrayElements(jb2, jpb2, 0);
	int len3 = env->GetArrayLength(jb3);
	jlong *jpb3 = env->GetLongArrayElements(jb3, 0);
	memcpy(b + 2 * pt.nWord, jpb3, len3 * sizeof(unsigned long long));
	env->ReleaseLongArrayElements(jb3, jpb3, 0);
	return pt.compute(b);
}

JNIEXPORT void JNICALL Java_phylonet_coalescent_Polytree_00024PTNative_cppBatchCompute
		(JNIEnv *env, jclass, jlongArray jres, jobjectArray jarr1, jobjectArray jarr2, jobjectArray jarr3){
	unsigned long long result[BATCH_SIZE] = {};
	unsigned long long* b = new unsigned long long[3 * BATCH_SIZE * pt.nWord]{};
	int size = env->GetArrayLength(jres);
	jlong *jpres = env->GetLongArrayElements(jres, 0);
	for (int i = 0; i < size; i += BATCH_SIZE){
		for (int j = 0; j < BATCH_SIZE && j < size - i; j++){
			jlongArray jb1 = (jlongArray) env->GetObjectArrayElement(jarr1, i + j);
			jlongArray jb2 = (jlongArray) env->GetObjectArrayElement(jarr2, i + j);
			jlongArray jb3 = (jlongArray) env->GetObjectArrayElement(jarr3, i + j);
			int len1 = env->GetArrayLength(jb1);
			jlong *jpb1 = env->GetLongArrayElements(jb1, 0);
			memcpy(b + (0 * BATCH_SIZE + j) * pt.nWord, jpb1, len1 * sizeof(unsigned long long));
			env->ReleaseLongArrayElements(jb1, jpb1, 0);
			int len2 = env->GetArrayLength(jb2);
			jlong *jpb2 = env->GetLongArrayElements(jb2, 0);
			memcpy(b + (1 * BATCH_SIZE + j) * pt.nWord, jpb2, len2 * sizeof(unsigned long long));
			env->ReleaseLongArrayElements(jb2, jpb2, 0);
			int len3 = env->GetArrayLength(jb3);
			jlong *jpb3 = env->GetLongArrayElements(jb3, 0);
			memcpy(b + (2 * BATCH_SIZE + j) * pt.nWord, jpb3, len3 * sizeof(unsigned long long));
			env->ReleaseLongArrayElements(jb3, jpb3, 0);
		}
		pt.compute(result, b);
		for (int j = 0; j < BATCH_SIZE && j < size - i; j++){
			jpres[i + j] = result[j];
		}
	}
	env->ReleaseLongArrayElements(jres, jpres, 0);
	delete b;
}

