package phylonet.coalescent;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;

import phylonet.tree.model.sti.STITreeCluster;

// TODO: why extend the abstract? It doesn't seem to follow the same pattern exactly
public class BipartitionWeightCalculator {

	AbstractInference inference;
	private AbstractDataCollection dataCollection;
	private int[] geneTreesAsInts;

	public BipartitionWeightCalculator(AbstractInference inference,
			int[] geneAsInts) {
		this.dataCollection = inference.dataCollection;
		this.inference = inference;
		this.geneTreesAsInts = geneAsInts;
	}

	
	class Intersects {
		long s0;
		long s1;
		long s2;
		long s3;


		public Intersects(long s0, long s1, long s2, long s3) {
			this.s0 = s0;
			this.s1 = s1;
			this.s2 = s2;
			this.s3 = s3;
		}
		
        public Intersects(Intersects side1, Intersects side2) {
            this(side1.s0+side2.s0,
					side1.s1+side2.s1,
					side1.s2+side2.s2,
					side1.s3+side2.s3);               
        }	

		public Intersects(Intersects other) {
			this.s0 = other.s0;
			this.s1 = other.s1;
			this.s2 = other.s2;
			this.s3 = other.s3;
		}

		public void addin(Intersects pop) {
			this.s0 += pop.s0;
			this.s1 += pop.s1;
			this.s2 += pop.s2;  
			this.s3 += pop.s3;
		}

		public void subtract(Intersects pop) {
			this.s0 -= pop.s0;
			this.s1 -= pop.s1;
			this.s2 -= pop.s2;         
			this.s3 -= pop.s3;   
		}

		public String toString() {
			return this.s0+","+this.s1+"|"+this.s2+","+this.s3;
		}
		
        public boolean isNotEmpty() {
        	return (this.s0 + this.s1 + this.s2 + this.s3) != 0;
        }
        
        public boolean hasEmpty() {
        	return this.maxPossible() == 0;
        }
        
        public long maxPossible() {
        	return (this.s0 * this.s1 * this.s2 * this.s3);
        }
	}

	private long allcases(Intersects side1, Intersects side2, Intersects side3) {
		return F(side1.s0,side2.s1,side3.s2,side3.s3)+
				F(side1.s1,side2.s0,side3.s2,side3.s3)+
				F(side1.s2,side2.s3,side3.s0,side3.s1)+
				F(side1.s3,side2.s2,side3.s0,side3.s1)+
				F(side3.s0,side2.s1,side1.s2,side1.s3)+
				F(side3.s1,side2.s0,side1.s2,side1.s3)+
				F(side3.s2,side2.s3,side1.s0,side1.s1)+
				F(side3.s3,side2.s2,side1.s0,side1.s1)+
				F(side1.s0,side3.s1,side2.s2,side2.s3)+
				F(side1.s1,side3.s0,side2.s2,side2.s3)+
				F(side1.s2,side3.s3,side2.s0,side2.s1)+
				F(side1.s3,side3.s2,side2.s0,side2.s1);
	}

	Intersects getSide(int i, Quadrapartition quart) {
		if (quart.cluster1.getBitSet().get(i)) {
			return new Intersects(1,0,0,0);
		} else if (quart.cluster2.getBitSet().get(i)) {
			return new Intersects(0,1,0,0);
		} else if (quart.cluster3.getBitSet().get(i)) {
			return  new Intersects(0,0,1,0);
		} else if (quart.cluster4.getBitSet().get(i)) {
			return  new Intersects(0,0,0,1);
		}
		else {
			return  new Intersects(0,0,0,0);
		}
	}
	
	public class Results {
		private double [] qs;
		private int effn;
		
		Results (double [] q, int n){
			setQs(q);
			setEffn(n);
		}

		public double [] getQs() {
			return qs;
		}

		public void setQs(double [] qs) {
			this.qs = qs;
		}

		public int getEffn() {
			return effn;
		}

		public void setEffn(int effn) {
			this.effn = effn;
		}
	}

	public Results getWeight(Quadrapartition [] quad ) {
		long [] fi = {0l,0l,0l};
		long mi = 0l;
		double [] weight = {0l,0l,0l};
		int effectiven = 0;
		Intersects []  allsides = new Intersects[3];
		boolean newTree = true;
		boolean cruise = false;
		
		//Quadrapartition []  quad = new Quadrapartition [] {quadm, new Quadrapartition(quadm.cluster1, quadm.cluster3, quadm.cluster2, quadm.cluster4), new Quadrapartition(quadm.cluster1, quadm.cluster4, quadm.cluster2, quadm.cluster3)};
		Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
		Deque<Intersects> [] stack = new Deque [] {new ArrayDeque<Intersects>(), new ArrayDeque<Intersects>(), new ArrayDeque<Intersects>()};

		for (Integer gtb: geneTreesAsInts){
			//n++;
			if (newTree) {
				STITreeCluster all = tit.next();
				for (int i=0; i<3; i++){
					allsides[i] = new Intersects(
						quad[i].cluster1.getBitSet().intersectionSize(all.getBitSet()),
						quad[i].cluster2.getBitSet().intersectionSize(all.getBitSet()),
						quad[i].cluster3.getBitSet().intersectionSize(all.getBitSet()),
						quad[i].cluster4.getBitSet().intersectionSize(all.getBitSet())
						);
				}
				newTree = false;
				mi = allsides[0].maxPossible();
				
				if ( mi != 0) {
					effectiven++;
				} else {
					cruise = true;
				}
				//sum +=  F(allsides.s0, allsides.s1, allsides.s2, allsides.s3);
			}
			if (gtb == Integer.MIN_VALUE) {
				if (!cruise) {
					//long fiall = fi[0] + fi[1] + fi[2];
					//if (fiall != 0)
					//double tf = 0.0; 
					for (int i=0; i<3; i++) {
						double efffreq = (fi[i]+0.0)/(2.0*mi);
						/*if (efffreq != 1 && efffreq != 0)
							System.err.println(efffreq);*/
						weight[i] += efffreq;
						//tf += efffreq;
					}
					//if (tf < .99) {
						//System.err.println("Warning: a gene tree is contributing only "+tf +" to total score of the branch with mi: "+mi);
					//}
				}
				for (int i=0; i<3; i++) stack[i].clear();
				newTree = true;
				cruise = false;
				fi = new long [] {0l,0l,0l};
				mi = 0;
			} else {
				if (cruise) continue;
				if (gtb >= 0){
					for (int i=0; i<3; i++) stack[i].push(getSide(gtb, quad[i]));
				} else  if (gtb == -2) {
					for (int i=0; i<3; i++) {
						Intersects side1 = stack[i].pop();
						Intersects side2 = stack[i].pop();
						Intersects newSide = new Intersects(side1, side2);
						stack[i].push(newSide);
						Intersects side3 = new Intersects(allsides[i]);
						side3.subtract(newSide);
		
						fi[i] += allcases(side1, side2, side3);
					}
	
					//geneTreesAsIntersects[n] = newSide;
				} else {
					for (int i=0; i<3; i++) {
						ArrayList<Intersects> children = new ArrayList<Intersects>();
					    Intersects newSide = new Intersects(0,0,0,0);
					    for (int j = gtb; j < 0 ; j++) {
					        Intersects pop = stack[i].pop();
					        children.add(pop);
					        newSide.addin(pop);
					    }
					    stack[i].push(newSide);
		                Intersects sideRemaining = new Intersects (allsides[i]);
		                sideRemaining.subtract(newSide);
		                if ( sideRemaining.isNotEmpty()) {
		                    children.add(sideRemaining);
		                }
		                
		                long sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0, sum01 = 0, sum23 = 0;
		                for (int j = 0; j < children.size(); j++){
		                	Intersects side = children.get(j);
		                	sum0 += side.s0;
		                	sum1 += side.s1;
		                	sum2 += side.s2;
		                	sum3 += side.s3;
		                	sum01 += side.s0 * side.s1;
		                	sum23 += side.s2 * side.s3;
		                }
		                for (int j = 0; j < children.size(); j++) {
		                	Intersects side = children.get(j);
		                	fi[i] += side.s0 * side.s1 * ((sum2 - side.s2) * (sum3 - side.s3) - sum23 + side.s2 * side.s3);
		                	fi[i] += side.s2 * side.s3 * ((sum0 - side.s0) * (sum1 - side.s1) - sum01 + side.s0 * side.s1);
		                }
					}
				}
			}
		}

		return  new Results(weight,effectiven);
	}
	
/*	private boolean checkFutileCalcs(Intersects side1, Intersects side2) {
		return ((side1.s0+side2.s0 == 0? 1 :0) +
    			(side1.s1+side2.s1 == 0? 1 :0) + 
    			(side1.s2+side2.s2 == 0? 1:0) > 1);
	}*/
	

	private long F(long a,long b,long c, long d) {
		if (a<0 || b<0 || c<0|| d<0) {
			throw new RuntimeException("negative side not expected: "+a+" "+b+" "+c);
		}
		return a*b*c*d;
	}	

	public class Quadrapartition {

		STITreeCluster cluster1;
		STITreeCluster cluster2;	
		STITreeCluster cluster3;
		STITreeCluster cluster4;
		private int _hash = 0;


		public Quadrapartition(STITreeCluster c1, STITreeCluster c2, STITreeCluster c3,STITreeCluster c4) {

			initialize(c1, c2, c3, c4);
		}
		
		private void initialize(STITreeCluster c1, STITreeCluster c2,
                STITreeCluster c3, STITreeCluster c4) {
            if (c1 == null || c2 == null || c3 == null) {
                throw new RuntimeException("none cluster" +c1+" "+c2+" "+c3);
            }
            int n1 = c1.getBitSet().nextSetBit(0), n2 = c2.getBitSet().nextSetBit(0), 
                    n3 = c3.getBitSet().nextSetBit(0), n4=c4.getBitSet().nextSetBit(0);
            int ntg1;
            int ntg2;
            STITreeCluster cluster_tmp1;
            STITreeCluster cluster_tmp2;    
            STITreeCluster cluster_tmp3;
            STITreeCluster cluster_tmp4;
            if (n1 < n2 ) {
                ntg1 = n1;
                cluster_tmp1 = c1;
                cluster_tmp2 = c2;
            }
            else {
                ntg1 = n2;
                cluster_tmp1 = c2;
                cluster_tmp2 = c1;
            }
            if (n3<n4) {
                ntg2 = n3;
                cluster_tmp3 = c3;
                cluster_tmp4 = c4;
            }
            else {
                ntg2 = n4;
                cluster_tmp3 = c4;
                cluster_tmp4 = c3;
            }
            
            if(ntg1<ntg2){
                cluster1 = cluster_tmp1;
                cluster2 = cluster_tmp2;
                cluster3 = cluster_tmp3;
                cluster4 = cluster_tmp4;
            }
            else{
                cluster1 = cluster_tmp3;
                cluster2 = cluster_tmp4;
                cluster3 = cluster_tmp1;
                cluster4 = cluster_tmp2;
            }
        }

		@Override
		public boolean equals(Object obj) {
			Quadrapartition trip = (Quadrapartition) obj; 

			return this == obj ||
					((trip.cluster1.equals(this.cluster1) && 
							trip.cluster2.equals(this.cluster2) && 
							trip.cluster4.equals(this.cluster4) && 
							trip.cluster3.equals(this.cluster3)));					
		}
		@Override
		public int hashCode() {
			if (_hash == 0) {
				_hash = cluster1.hashCode() * cluster2.hashCode()
						* cluster4.hashCode() * cluster3.hashCode();
			}
			return _hash;
		}
		@Override
		public String toString() {		
			return cluster1.getBitSet().toString2()+"|"+cluster2.getBitSet().toString2()+
					"#"+cluster3.getBitSet().toString2()+"|"+cluster4.getBitSet().toString2();
		}
		public String toString2() {
			return cluster1.toString()+"|"+cluster2.toString()+
					"#"+cluster3.toString()+ "|"+cluster4.toString();
		}


	}

}
