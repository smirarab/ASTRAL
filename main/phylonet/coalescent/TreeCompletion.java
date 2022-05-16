package phylonet.coalescent;

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import com.sun.org.apache.bcel.internal.generic.NEW;

import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.TMutableNode;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.util.Trees;

class NodeInfo{
	private String color="";	
	private int greenDescendents = 0;
	public NodeInfo(String color, int greenDes){
		this.color = color;
		this.greenDescendents = greenDes;
	}
	public NodeInfo(int greenDes){
		this.greenDescendents = greenDes;
	}
	
	public String getColor() {
		return color;
	}
	public int getGrDes() {
		return greenDescendents;
	}
	
}
public class TreeCompletion {
	
	public static void main(String[] args) throws IOException, ParseException{
		//String tr1 = "((1,(2,3)a)b,(4,(5,(6,7)d)c));";
//		String tr1 = "((1,2)a,(3,(4,(6,(5,7))d)c));";
		//((3,99),((11,12,(22,1),30),(4,((5,9,18)a,6)))b,(7,2))c;

//		String tr1 = "((7,2)a,(3,(4,(6,(5,(9,(18,23,90))))d)c)f,(11,12,(22,1),30)g,99)r;";
//		String tr2 = "((3,99),((11,22,1),(5,9,18)a)b)c;";
		
//		String tr2 = "(193,(15,(62,(((99,43),(26,(106,(137,78)))),(((66,(3,(102,(151,165)))),(16,((186,((76,(124,(51,79))),(98,(123,((107,195),(131,(168,104))))))),(((95,121),(9,136)),((6,(67,41)),(81,((45,130),(198,(12,47))))))))),((((153,53),(159,44)),((183,42),((196,22),((19,(40,(14,110))),(156,((173,111),(13,185))))))),((83,((74,148),(146,(109,118)))),(((8,(50,10)),((24,88),((92,38),(132,(75,154))))),((31,((134,(157,82)),((49,72),((188,149),(119,(184,143)))))),(128,(((69,(((29,158),(71,96)),((179,28),((100,87),((126,122),(59,(120,23))))))),((181,((77,30),(152,(70,(85,46))))),((150,(187,(37,61))),((89,182),(161,(135,7)))))),(0,((((113,36),(63,34)),((172,178),((39,65),(68,(166,(60,86)))))),(((145,48),((127,(169,(171,(199,176)))),((32,(18,191)),(180,((142,(105,(117,175))),(((17,93),((162,52),(108,54))),((20,94),((138,(192,160)),(197,(200,141)))))))))),((((58,(114,167)),(57,(139,25))),(2,(((33,1),(190,21)),(5,((80,163),(55,(101,(189,84)))))))),(((73,112),(56,(4,140))),(91,((147,11),(((97,133),((177,27),(144,129))),((64,(125,((90,170),(174,(164,103))))),(35,((115,194),(155,116)))))))))))))))))))))));";
//		String tr1 = "(193,(15,(62,(((99,43),(26,(106,(137,78)))),(((66,(3,(102,(151,165)))),(16,((186,((76,(124,(51,79))),(98,(123,((107,195),(131,(168,104))))))),(((95,121),(9,136)),((6,(67,41)),(81,((45,130),(198,(12,47))))))))),((((153,53),(159,44)),((183,42),((196,22),((19,(40,(14,110))),(156,((173,111),(13,185))))))),((83,((74,148),(146,(109,118)))),(((8,(50,10)),((24,88),((92,38),(132,(75,154))))),((31,((134,(157,82)),((49,72),((188,149),(119,(184,143)))))),(128,(((69,(((29,158),(71,96)),((179,28),((100,87),((126,122),(59,(120,23))))))),((181,((77,30),(152,(70,(85,46))))),((150,(187,(37,61))),((89,182),(161,(135,7)))))),(0,((((113,36),(63,34)),((172,178),((39,65),(68,(166,(60,86)))))),(((145,48),((127,(169,(171,(199,176)))),((32,(18,191)),(180,((142,(105,(117,175))),(((17,93),((162,52),(108,54))),((20,94),((138,(192,160)),(197,(200,141)))))))))),((((58,(114,167)),(57,(139,25))),(2,(((33,1),(190,21)),(5,((80,163),(55,(101,(189,84)))))))),(((73,112),(56,(4,140))),(91,((147,11),(((97,133),((177,27),(144,129))),((64,(125,((90,170),(174,(164,103))))),(35,((115,194),(155,116)))))))))))))))))))))));";
		//(18,(7,2),((3,99),4),((11,(5,9)a,0,22,30)b,6))c;
//		String tr2 = "(9,7,3,4);";
		String tr1 = args[0];
		String tr2 = args[1]; //backbone
		NewickReader nr = new NewickReader(new StringReader(tr1));
		STITree<Double> gt = new STITree<Double>(true);
		nr.readTree(gt);
		
		NewickReader nr2 = new NewickReader(new StringReader(tr2));
		STITree<Double> st = new STITree<Double>(true);
		nr2.readTree(st);
		STINode newNode = new STITree(((STINode<Double>) gt.getNode("18")).toNewick()).getRoot();

		System.err.println(st.toNewick());
		st = treeCompletion(gt,st);
		System.err.println("out: "+st.toNewick());

	}
	
//	static STITree addToTree(STITree tree , STINode adoptingNode, STINode toMoveNode){
//		
//		STINode newNode = null;
//		try {
//			newNode = new STITree(((STINode<Double>) toMoveNode).toNewick()).getRoot();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (ParseException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		if(!adoptingNode.isRoot()){
//			STINode newinternalnode = adoptingNode.getParent().createChild();
//			newinternalnode.adoptChild(adoptingNode);
//			newinternalnode.createChild(newNode);
//			return tree;
//		}
//		else{
//
//			STINode newinternalnode = adoptingNode.createChild();
//			TNode child = (TNode) adoptingNode.getChildren().iterator().next();
//			newinternalnode.adoptChild((TMutableNode) child );
//			newinternalnode.createChild(newNode);
//			return tree;
//			
//		}
//		
//	}

	
static STITree addToTreePolytomy2(STITree tree , STINode adoptingNode, ArrayList<STINode> redChildren){
		
		ArrayList<STINode> newnodes = new ArrayList<STINode>();
		try {
			for(STINode n : redChildren){
				if(!n.isLeaf())
					n.setName("");
				newnodes.add(new STITree(((STINode<Double>)n).toNewick()).getRoot());
		}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if(adoptingNode.isLeaf()){			
			for(STINode n : newnodes){
				adoptingNode.getParent().createChild(n);
			}
		}
		else{
			for(STINode n : newnodes){
				adoptingNode.createChild(n);
			}			
		}
		return tree;
	}

static STITree addToTreePolytomy(STITree tree , STINode adoptingNode, ArrayList<STINode> redChildren){
		
		ArrayList<STINode> newnodes = new ArrayList<STINode>();
		try {
			for(STINode n : redChildren){
				if(!n.isLeaf())
					n.setName("");
				newnodes.add(new STITree(((STINode<Double>)n).toNewick()).getRoot());
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		if(!adoptingNode.isRoot()){
			
			STINode newinternalnode = adoptingNode.getParent().createChild();
			newinternalnode.adoptChild(adoptingNode);
			for(STINode n : newnodes){
				newinternalnode.createChild(n);
			}			
			
		}
		else{
				
			STINode newinternalnode = adoptingNode.createChild();
			TNode child = (TNode) adoptingNode.getChildren().iterator().next();
			newinternalnode.adoptChild((TMutableNode) child );
			for(STINode n : newnodes)
				newinternalnode.createChild(n);
			
		}
		return tree;
		
	}
	
	static void nodeColoring(Tree stTree, Tree gtTree) {
		
		
		for (TNode node : stTree.postTraverse()) {
			if(node.isLeaf())
				((STINode) node).setData(new NodeInfo(1));
			else{
				int greenCount = 0;
				for (TNode child:node.getChildren())		
					greenCount += ((NodeInfo) ((STINode) child).getData()).getGrDes();
				((STINode) node).setData(new NodeInfo(greenCount));	
		
			}
		}
			
		Set<String> stLeaves = new HashSet<String>(Arrays.asList(stTree.getLeaves())); 
		
		for (TNode gtNode : gtTree.postTraverse()) {
			if(gtNode.isLeaf()){
				String name = gtNode.getName();
				if (stLeaves.contains(name)){
					((STINode) gtNode).setData(new NodeInfo("G", 1));
				}
				else{
					((STINode) gtNode).setData(new NodeInfo("R", 0));
				}
			}
			else{
				boolean allred = true;
				boolean allgreen = true;
				boolean hasred = false;
				int greenCount = 0;
				for (TNode child: gtNode.getChildren()){
					
					NodeInfo info = (NodeInfo) ((STINode) child).getData();
					String data = info.getColor();
					greenCount += info.getGrDes();
//					String data = (String) ((STINode) child).getData();
					if(data.equals("B") || data.equals("BM") ){
						allred = false;
						allgreen = false;
					}
					else if(data.equals("G")){
						allred = false;
					}
					else if(data.equals("R")){
						allgreen = false;
						hasred = true;
					}
				}
					if(allred) 
						((STINode) gtNode).setData(new NodeInfo("R", greenCount));
					else if(allgreen)
						((STINode) gtNode).setData(new NodeInfo("G", greenCount));
					else if(hasred)
						((STINode) gtNode).setData(new NodeInfo("BM", greenCount));
					else
						((STINode) gtNode).setData(new NodeInfo("B", greenCount));
			}
		}
		
	}
	
	static ArrayList<STITree> treeCompletionRepeat(STITree referenceTree, STITree constraintTree, int REPEATS){
		SchieberVishkinLCA lcaLookup = new SchieberVishkinLCA(constraintTree);
		String[] gtLeaves = referenceTree.getLeaves();
		String[] stLeaves = constraintTree.getLeaves();
		List<String> common = new ArrayList<String>(Arrays.asList(gtLeaves));
		common.retainAll(Arrays.asList(stLeaves));
		ArrayList<STITree> results = new ArrayList<STITree>();
		ArrayList<STITree> temps = new ArrayList<STITree>();
		for(int i=0; i< REPEATS; i++){
			temps.add(new STITree(constraintTree));
		}
		
		ArrayList<Integer> randomRoots = new ArrayList<Integer>();
		for(int i=0;i< REPEATS && i < common.size() ;i++){
			randomRoots.add(GlobalMaps.random.nextInt(common.size()));
		}

		for(int i=0 ; i< REPEATS && i < common.size(); i++){
			String root = common.get(randomRoots.get(i));
			//System.err.println("Both constraint tree and gene tree repeat "+i+" are rooted at "+root);
			Trees.removeBinaryNodes(temps.get(i));
			temps.get(i).rerootTreeAtEdge(temps.get(i).getNode(root));
			Trees.removeBinaryNodes(referenceTree);
			referenceTree.rerootTreeAtEdge(referenceTree.getNode(root));
			results.add(treeCompletion(referenceTree, temps.get(i)));
		}
		return results;
		
	}
	static STITree treeCompletion(STITree gTree, STITree sTree){
		nodeColoring(sTree, gTree);
		HashMap<Integer, Integer> LCAMap = createLCAMap(sTree, gTree);
		LCAMap = resolvePolytomies(sTree, gTree, LCAMap);
	  
	        // Create an empty stack and push root to it 
	        Stack<TNode> nodeStack = new Stack<TNode>(); 
	        nodeStack.push(gTree.getRoot());
	        while (nodeStack.empty() == false) { 
	              
	            // Pop the top item from stack and print it 
	            TNode mynode = nodeStack.peek();  
	            NodeInfo info = (NodeInfo) ((STINode) mynode).getData();
				String data = info.getColor();

	            if(data.equals("BM")){
	            	
//	            	int childrenCount = mynode.getChildCount();
	            	int redchild = 0;
	            	ArrayList<STINode> redChildren = new ArrayList<STINode>();
	            	ArrayList<STINode> nonRed = new ArrayList<STINode>();
	            	for(TNode child:mynode.getChildren()){
						String childData = (String) ((NodeInfo) ((STINode) child).getData()).getColor();
	            		if(childData.equals("R")){
	            			redChildren.add((STINode) child);            			
	            			redchild += 1;
	            		}
	            		else{
	            			nonRed.add((STINode) child);
	            		}
	            	}
	            	int id = LCAMap.get(mynode.getID());
	    	        STINode snode = sTree.getNode(id);
	    	        
	    	        if(nonRed.size()==1 && ((NodeInfo) snode.getData()).getGrDes() == ((NodeInfo) nonRed.get(0).getData()).getGrDes()){
	    	        		sTree = addToTreePolytomy(sTree, snode,  redChildren); 
	    	        }
	    	        else{
	    	        		sTree = addToTreePolytomy2(sTree, snode,  redChildren);// does not create a new node
	    	        }
	    	        
	    	        //if there is only one non-red child,it means that we definitely don't have this 
	    	        //node in the other tree and we should create a new node
////	    	        if(redchild >= 1 && mynode.getChildCount()-redchild == 1){	
////	    	        	sTree = addToTreePolytomy(sTree, snode,  redChildren);
////	    	        }
//	    	        if(redchild ==1 && mynode.getChildCount() == 2){	            	
//		    	        sTree = addToTree(sTree, snode, (STINode) redChildren.get(0));		    			
//	            	}
//	            	else{
//	            		sTree = addToTreePolytomy(sTree, snode,  redChildren);
//	            	}
	    	        
	    	        //----------------------------
				}
	            nodeStack.pop(); 
	            
	            // Push children of current node to the stack
	            //since at first we added children of the root, we skip it here

	            for(TNode t: mynode.getChildren())
	            		nodeStack.push(t); 

	    }  
	        return sTree;
	}
	
	
	static HashMap<Integer,Integer> resolvePolytomies(Tree stTree, Tree gtTree, HashMap<Integer,Integer> LCAMap){
		for (TNode node : gtTree.postTraverse()) {
			if (node.isLeaf()) {
			}
			else{
            	ArrayList<STINode> nonRed = new ArrayList<STINode>();
            	for(TNode child: node.getChildren()){
					String childData = (String) ((NodeInfo) ((STINode) child).getData()).getColor();
            		if(!childData.equals("R"))     			
            			nonRed.add((STINode) child);
            		
            	}
            	if(nonRed.size() >= 2){
            		STINode snode = (STINode) stTree.getNode(LCAMap.get(node.getID()));
            		//need to add resolution
            		if(snode.getChildCount() > nonRed.size()){
            			boolean sameLCAs = true;
            			ArrayList<STINode> LCAchildren = new ArrayList<STINode>();
            			for(TNode child: nonRed){
            				STINode snodechild = (STINode) stTree.getNode(LCAMap.get(child.getID()));
            				if(!snodechild.getParent().equals(snode)){
            					sameLCAs = false;
            				}
            				else{
            					LCAchildren.add(snodechild);
            				}
            			}
            			if(sameLCAs){
            				STINode newinternalnode = snode.createChild();
            				for(TNode child: LCAchildren){
            					newinternalnode.adoptChild((TMutableNode) child);
            				}
            				LCAMap.put(node.getID(), newinternalnode.getID());
            			}
            			
            		}
            	}
			}
		}
		return LCAMap;
	}
	
	static HashMap<Integer,Integer> createLCAMap(Tree stTree, Tree gtTree) {

		HashMap<Integer, Integer> LCAMap = new HashMap<Integer, Integer>();
		SchieberVishkinLCA lcaLookup = new SchieberVishkinLCA(stTree);

		Stack<TNode> stack = new Stack<TNode>();
		for (TNode gtNode : gtTree.postTraverse()) {
					if (gtNode.isLeaf()) {

						String color = (String) ((NodeInfo) ((STINode) gtNode).getData()).getColor();
						
						if(color.equals("G")){
							TNode t = stTree.getNode(gtNode.getName());
							stack.push(t);
//							LCAMap.put(gtNode.getID(), new ArrayList<Integer>(Arrays.asList(t.getID(),0)));
							LCAMap.put(gtNode.getID(), t.getID());
						}
						if(color.equals("R")){
							stack.push(null);
						}
					} else {
						if(gtNode.getChildCount()==2){
							TNode rightLCA = stack.pop();
							TNode leftLCA = stack.pop();
							// If gene trees are incomplete, we can have this case
							if (rightLCA == null && leftLCA == null) {
								stack.push(null);
								continue;
							}
							else if (rightLCA == null || leftLCA == null) {
								if(rightLCA != null){
									stack.push(rightLCA);
									LCAMap.put(gtNode.getID(), rightLCA.getID());
								}
								if(leftLCA != null){
									stack.push(leftLCA);
									LCAMap.put(gtNode.getID(), leftLCA.getID());
								}
	
								continue;
							}
							
							TNode lca = lcaLookup.getLCA(leftLCA, rightLCA);
							stack.push(lca);
							LCAMap.put(gtNode.getID(), lca.getID());
						}
						else{
							TNode[] children = new TNode[gtNode.getChildCount()];
							boolean allnull = true;
							for(int i =0;i < gtNode.getChildCount();i++){
								children[i] = stack.pop();
								if(children[i] != null)
									allnull = false;
							}
							if (allnull) {
								stack.push(null);
								continue;
							}
							else {
								Set<TNode> notnulls = new HashSet<TNode>();
								for(TNode child: children){
									if(child != null)
										notnulls.add(child);
								}
								TNode lca = lcaLookup.getLCA(notnulls);								
								stack.push(lca);
								LCAMap.put(gtNode.getID(), lca.getID());
								continue;
							}
							
						}

					}
				}
				return LCAMap;
			}
	
	}


	

