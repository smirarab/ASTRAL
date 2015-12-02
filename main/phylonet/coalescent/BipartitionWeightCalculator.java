package phylonet.coalescent;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;

class BipartitionWeightCalculator extends AbstractWeightCalculator<Tripartition> {

	WQInference inference;
	private WQDataCollection dataCollection;

	public BipartitionWeightCalculator(AbstractInference<Tripartition> inference) {
		super(false);
		this.dataCollection = (WQDataCollection) inference.dataCollection;
		this.inference = (WQInference) inference;
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

	public long getWeight(Quadrapartition quad) {
		long weight = 0l;
		Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();

		Deque<Intersects> stack = new ArrayDeque<Intersects>();

		Intersects  allsides = null;
		boolean newTree = true;

		for (Integer gtb: dataCollection.geneTreesAsInts){
			//n++;
			if (newTree) {
				STITreeCluster all = tit.next();
				allsides = new Intersects(
						quad.cluster1.getBitSet().intersectionSize(all.getBitSet()),
						quad.cluster2.getBitSet().intersectionSize(all.getBitSet()),
						quad.cluster3.getBitSet().intersectionSize(all.getBitSet()),
						quad.cluster4.getBitSet().intersectionSize(all.getBitSet())
						);
				newTree = false;
				//sum +=  F(allsides.s0, allsides.s1, allsides.s2, allsides.s3);
			}
			if (gtb >= 0){
				stack.push(getSide(gtb, quad));
			} else if (gtb == Integer.MIN_VALUE) {
				stack.clear();
				newTree = true;
			} else if (gtb == -2) {
				Intersects side1 = stack.pop();
				Intersects side2 = stack.pop();
				Intersects newSide = new Intersects(side1, side2);
				stack.push(newSide);
				Intersects side3 = new Intersects(allsides);
				side3.subtract(newSide);

				weight+= allcases(side1, side2, side3);

				//geneTreesAsIntersects[n] = newSide;
			} else {
				/*Intersects side1 = stack.pop();
				Intersects side2 = stack.pop();
				Intersects side3 = stack.pop();
				Intersects newSide = new Intersects(
						side1.s0+side2.s0+side3.s0,
						side1.s1+side2.s1+side3.s1,
						side1.s2+side2.s2+side3.s2,
						side1.s3+side2.s3+side3.s3);
				stack.push(newSide);
				Intersects rem = new Intersects(allsides);
				rem.subtract(newSide);

				if (rem.s0+rem.s1+rem.s2+rem.s3 >0) {
					throw new RuntimeException("polytomies are too complicated :( Ask us later");						
				}
				weight+= allcases(side1, side2, side3);*/
				ArrayList<Intersects> children = new ArrayList<Intersects>();
			    Intersects newSide = new Intersects(0,0,0,0);
			    for (int i = gtb; i < 0 ; i++) {
			        Intersects pop = stack.pop();
			        children.add(pop);
			        newSide.addin(pop);
			    }
			    stack.push(newSide);
                Intersects sideRemaining = new Intersects (allsides);
                sideRemaining.subtract(newSide);
                if ( sideRemaining.isNotEmpty()) {
                    children.add(sideRemaining);
                }
                for (int i = 0; i < children.size(); i++) {
                    Intersects side1 = children.get(i);
                    
                    for (int j = i+1; j < children.size(); j++) {
                        Intersects side2 = children.get(j);
                        //if (children.size() > 5 && checkFutileCalcs(side1,side2))
                        //	continue;
                        
                        for (int k = j+1; k < children.size(); k++) {
                            Intersects side3 = children.get(k);
                            weight += allcases(side1,side2,side3);
                        }
                    }
                }
			}
		}

		return  weight/2;
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

	class Quadrapartition {

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
			cluster1 = c1;
			cluster2 = c2;	
			cluster3 = c3;
			cluster4 = c4;
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
			return cluster1.toString()+"|"+cluster2.toString()+
					"||"+cluster3.toString()+cluster4.toString();
		}


	}

	@Override
	protected void prepareWeightTask(
			ICalculateWeightTask<Tripartition> weigthWork,
			AbstractComputeMinCostTask<Tripartition> task) {
		// TODO Auto-generated method stub

	}

	@Override
	public ICalculateWeightTask<Tripartition> getWeightCalculateTask(
			Tripartition t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees) {
		// TODO Auto-generated method stub

	}

}