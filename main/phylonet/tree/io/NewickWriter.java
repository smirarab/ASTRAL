/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package phylonet.tree.io;

import java.io.PrintWriter;
import java.io.Writer;
import java.util.Iterator;

import phylonet.tree.model.TMutableNode;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;



/**
 * This class writes trees in newick format to a writer.
 * 
 * @author Derek Ruths
 */
public class NewickWriter {

	// fields
	private PrintWriter _writer;
	
	// constructors
	public NewickWriter(Writer w) {
		_writer = new PrintWriter(w);
	}
	
	// methods
	public void writeTree(TNode node, boolean end_line) {
		writeNode(node);
		
		_writer.print(";");
		
		if(end_line) {
			_writer.println();
		}
		
		_writer.flush();
	}
	
	public void writeTreeWD(TNode node, boolean end_line) {
		writeNodeWD(node);
		
		_writer.print(";");
		
		if(end_line) {
			_writer.println();
		}
		
		_writer.flush();
	}
	
	public void writeTree(Tree tree, boolean end_line) {
		
		if(tree.getRoot() != null) {
			writeNode(tree.getRoot());
		}
		
		_writer.print(";");
		
		if(end_line) {
			_writer.println();
		}
		
		_writer.flush();
	}
	
	public void writeTreeWD(Tree tree, boolean end_line) {
		
		if(tree.getRoot() != null) {
			writeNodeWD(tree.getRoot());
		}
		
		_writer.print(";");
		
		if(end_line) {
			_writer.println();
		}
		
		_writer.flush();
	}
	
	private void writeNode(TNode node) {
		
		if(node.getChildCount() > 0) {
			_writer.print("(");

			Iterator<? extends TNode> i = node.getChildren().iterator();
			while(i.hasNext()) {
				TNode child = i.next();
				
				writeNode(child);
				
				if(i.hasNext()) {
					_writer.print(",");
				}
			}
			
			_writer.print(")");
		}
		
		// write this node's info
		_writer.print(node.getName());
		if(node.getParentDistance() != TMutableNode.NO_DISTANCE) {
			_writer.print(":" + node.getParentDistance());
			
		}
	}
	
	private void writeNodeWD(TNode node) {
		
		if(node.getChildCount() > 0) {
			_writer.print("(");

			Iterator<? extends TNode> i = node.getChildren().iterator();
			while(i.hasNext()) {
				TNode child = i.next();
				
				writeNodeWD(child);
				
				if(i.hasNext()) {
					_writer.print(",");
				}
			}
			
			_writer.print(")");
		}
		
		// write this node's info	
		_writer.print(node.getName());
		if(((STINode)node).getData() != null) {
			_writer.print(((STINode)node).getData().toString());
		}
		if(node.getParentDistance() != TMutableNode.NO_DISTANCE) {
			_writer.print(":" + node.getParentDistance());
			
		}
		
	}
}
