package phylonet.coalescent;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;


public class SimilarityMatrixMP extends SimilarityMatrix{

	public SimilarityMatrixMP(float[][] from) {
		super(from);
	}

	public SimilarityMatrixMP(int n) {
		super(n);
	}

	@Override
	protected int getMemChunkCount() {
		return Threading.getDistMatrixChunkSize();
	}

	@Override
	protected int getChunkSize(final List<Tree> geneTrees) {
		return (int) Math.ceil(geneTrees.size()/(Threading.getNumThreads()+0.0));
	}

	Object[][][] locks = null;


	@Override
	protected void updateDnArray(int m, Long[][] dn, int c, int l, int r) {
		int ll = l <=r ? l : r;
		int rr = l <=r ? r : l;
		synchronized (locks[m][ll][rr]) {
			dn[l][r] += c*(c-1)/2;
			dn[r][l] = dn[l][r];
		}
	}

	@Override
	protected void updateArray(int m, Long[][] array, long sim, int l, int r) {
		int ll = l <=r ? l : r;
		int rr = l <=r ? r : l;
		synchronized (locks[m][ll][rr]) {
			array[l][r] += sim;
			array[r][l] = array[l][r];
		}
	}

	@Override
	public void populateByQuartetDistance(final List<STITreeCluster> treeAllClusters,final List<Tree> geneTrees) {
		Logging.log("with " + getMemChunkCount() + " distance matrices for parallellism");

		this.matrix = new float[n][n];
		Logging.logTimeMessage("SimilarityMatrix 145-148: ");

		final int chunksize = getChunkSize(geneTrees);
		final int memchunkcount = getMemChunkCount();
		//final int memchunksize = (int) Math.ceil((geneTrees.size()+0.0)/memchunkcount);

		final Long[][][] a = new Long[memchunkcount][n][n];
		final Long[][][] d = new Long[memchunkcount][n][n];
		locks = new Object[getMemChunkCount()][n][n];
		for (int c = 0 ; c < memchunkcount; c++) {
			for(int i = 0; i < n; i++) {
				for(int j = 0; j < n; j++) {
					a[c][i][j]  = 0l;
					d[c][i][j] = 0l;
					locks[c][i][j] = new Object();
				}
			}
		}


		ArrayList<Future> futures = new ArrayList<Future>();
		for (int i = 0; i < geneTrees.size()/chunksize + 1; i+= 1) {
			final int consti = i;
			futures.add(Threading.submit( new Callable<Boolean>() {
				public Boolean call() {
					int start = consti * chunksize;
					int end = Math.min(start + chunksize, geneTrees.size());
					int m = consti % memchunkcount;	
					Long[][] array = a[m];
					Long[][] dn = d[m];
					for(int w = start; w < end; w++) {
						processGene(m, array, dn, treeAllClusters.get(w), geneTrees.get(w));
					}
					return Boolean.TRUE;
				}
			}

					));
		}
		for (Future future: futures) {
			try {
				future.get();
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}
		locks = null;
		for (int c = 0 ; c < memchunkcount; c++) {
			for(int i = 0; i < n; i++) {
				for(int j = 0; j < n; j++) {
					a[0][i][j] += a[c][i][j];
					//System.err.print(a[c][i][j]);
					//System.err.print(",");
					d[0][i][j] += d[c][i][j];
				}
				//Logging.log();
			}
			//Logging.log();
		}
		Logging.logTimeMessage("SimilarityMatrix 161-164: ");

		normalize(a[0],d[0]);

		/*Logging.log();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				System.err.print(String.format("%.5f",matrix[i][j]));
				System.err.print(",");
			}
			Logging.log();
		}*/

	}
}
