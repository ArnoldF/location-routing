package MDVRP;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;


public class LinKernighan {

	private static class NOptMove {
		public ArrayList<Edge> brokenEdges;
		public ArrayList<Edge> addedEdges;
		public double saving;
		public double oldLength;
		public int nextNode;
		public int degree[];
		
		public NOptMove()
		{
			brokenEdges = new ArrayList<Edge>();
			addedEdges = new ArrayList<Edge>();
			saving = 0d;
			oldLength = 0d;
		}
		public NOptMove(double s, double o, int n, ArrayList<Edge> b, ArrayList<Edge> a, int[] d)
		{
			brokenEdges = b;
			addedEdges = a;
			saving = s;
			oldLength = o;
			nextNode = n;
			degree = d;
		}
		
		public NOptMove copy()
		{
			ArrayList<Edge> b = new ArrayList<Edge>(); 
			ArrayList<Edge> a = new ArrayList<Edge>();
			for (int i=0;i<addedEdges.size();i++)
				a.add( addedEdges.get(i).copy() );
			for (int i=0;i<brokenEdges.size();i++)
				b.add( brokenEdges.get(i).copy() );
			return new NOptMove(saving,oldLength,nextNode,b,a,null);	
		}
	}

	public Route linKernighan(Route r, int lambda, MDVRPModel model)
	{
		//-------------------------------------------------------------------
		// INITIALISATION
		//-------------------------------------------------------------------
		int num_nodes = r.length+1;
		NOptMove bestMove = new NOptMove();
		boolean improved = true;
		int startIndex = 0;
		int count = 0;
		int routeIndex = r.nodes.get(1).routeIndex;
		
		//for each node, create a list of candidate neighbors (10 nearest, sorted by distance)
		HashMap<Integer,ArrayList<Edge>> candidate_neighbours = new HashMap<Integer,ArrayList<Edge>>();
		for (int i=0; i<num_nodes; i++)
		{
			ArrayList<Edge> distance = new ArrayList<Edge>();
			for (int j=0; j<num_nodes; j++)
				if (i!=j)// && (r.nodes.get(i).customers_considered_as_neigbours.contains( r.nodes.get(j).index ) || r.nodes.get(i).index >= model.customers ) )
					distance.add( new Edge( model.getDistance(r.nodes.get(i).index, r.nodes.get(j).index),i,r.nodes.get(j).index ) );
			
			//get the 10 nearest nodes, sorted by distance
			Collections.sort( distance );
			Collections.reverse( distance );
			int N = Math.max(0, distance.size() - 10);
			distance.subList(distance.size() - N, distance.size()).clear();
			candidate_neighbours.put(r.nodes.get(i).index,distance);
		}
		
		//keep track of the position of each node within the route
		HashMap<Integer,Integer> position = new HashMap<Integer,Integer>();
		for (int i=0; i<num_nodes; i++)
			position.put(r.nodes.get(i).index, i);
		int[] nodeAt = new int[num_nodes];
		for (int i=0; i<num_nodes; i++)
			nodeAt[i] = r.nodes.get(i).index;
		
		//-------------------------------------------------------------------
		// EXECUTE AS LONG AS IMPROVING MOVES ARE FOUND
		//-------------------------------------------------------------------
		mainloop2 : 
		while ( improved && count < 1000)
		{			
			improved = false;
			count++;
			
			//order and safe all current edges
			ArrayList<Edge> current_edges = new ArrayList<Edge>();
			for (int i=0; i<num_nodes-1; i++)
				current_edges.add( new Edge( model.getDistance(r.nodes.get(i).index, r.nodes.get(i+1).index), i, i+1));
			current_edges.add( new Edge( model.getDistance(r.nodes.get(0).index, r.nodes.get(num_nodes-1).index), 0, num_nodes-1 ));
			Collections.sort( current_edges );

			//starting with the removal of each edge, try to build a move
			for (int startEdge=startIndex; startEdge<startIndex+current_edges.size(); startEdge++)
			//for (int dir=0;dir<2;dir++)
			{
				int index;
				if (startEdge >= current_edges.size())
					index = startEdge - current_edges.size();
				else
					index = startEdge;
				ArrayList<NOptMove> NOptList = new ArrayList<NOptMove>();
				NOptMove bestNOptMove = new NOptMove();
				int depth=1;	
					
				//keep track of how many incident edges each node has
				int[] degree = new int[num_nodes];
				for (int i=0; i<num_nodes; i++)
					degree[i] = 2;
				
				//---------------------------------------------
				//remove initial edge
				//---------------------------------------------
				int nextNode = current_edges.get(index).node1;
				int endNode	 = current_edges.get(index).node2;
				//if (dir==1)
				{
					//nextNode = current_edges.get(index).node2;
					//endNode	 = current_edges.get(index).node1;
				}
				ArrayList<Edge> brokenEdges = new ArrayList<Edge>();
				ArrayList<Edge> addedEdges  = new ArrayList<Edge>();				
				degree[nextNode]--;
				degree[endNode]--;
				brokenEdges.add(new Edge(nextNode, endNode));
				double oldLength = model.getDistance(r.nodes.get(nextNode).index, r.nodes.get(endNode).index);
				double newLength = 0d;
				double saving = 0d;
				NOptMove move = new NOptMove(saving, oldLength, nextNode, brokenEdges, addedEdges, degree);
				NOptList.add( move );
				
				mainloop:
				while (depth<=lambda)					
				{
					//---------------------------------------------
					//try to close all incompleted moves by adding an edge to endNode
					//---------------------------------------------
					if (depth>1)
					{
						boolean improvement = false;
						for (int m=0; m<NOptList.size(); m++)
						{
							nextNode 	= NOptList.get(m).nextNode;
							newLength 	= model.getDistance( r.nodes.get(nextNode).index, r.nodes.get(endNode).index );
							saving		= NOptList.get(m).saving;
							saving 		+= NOptList.get(m).oldLength - newLength;
							model.ev1++;
							
							//if the move improves the route, safe it as candidate move
							if ( saving > bestNOptMove.saving && !areNeighbours(nextNode, endNode, num_nodes) )
							{
								ArrayList<Edge> addedEdges_copy = copyEdgeList(NOptList.get(m).addedEdges);
								addedEdges_copy.add( new Edge( nextNode, endNode ) );
								if (!hasSubTour(addedEdges_copy, NOptList.get(m).brokenEdges, num_nodes) )
								//we found a better move
								{
									bestNOptMove.saving = saving;
									bestNOptMove.addedEdges = addedEdges_copy;
									bestNOptMove.brokenEdges = copyEdgeList(NOptList.get(m).brokenEdges);
									improvement = true;
								}
							}
						}
					}//closure
					
					//---------------------------------------------
					//continue all incomplete moves with positive partial gain
					//---------------------------------------------
					if (depth < lambda)
					{
						int num_partial_moves = NOptList.size();
						for (int m=0; m<num_partial_moves; m++)
						{
							//find a new edge
							nextNode = NOptList.get(0).nextNode;
							int consideredNeighbours = candidate_neighbours.get( nodeAt[nextNode]).size();
							//if (depth>=2)
								//consideredNeighbours = Math.min(2, sortedDistance.get(nextNode).size());
						
							//try to connect the current node with all candidate neighbours
						    for (int candidate=0; candidate<consideredNeighbours; candidate++)
						    if ( !areNeighbours(nextNode, position.get( candidate_neighbours.get( nodeAt[nextNode] ).get(candidate).node2), num_nodes) )
							{
								int new_neighbour = position.get( candidate_neighbours.get( nodeAt[nextNode] ).get(candidate).node2 );
								boolean neighbour_found = false;
								saving	 = NOptList.get(0).saving;
								degree   = NOptList.get(0).degree.clone();
								if ( new_neighbour != endNode && degree[new_neighbour]>0 && !containsEdge(NOptList.get(0).addedEdges,new Edge(nextNode,new_neighbour)) )
								{
									neighbour_found = true;
									newLength	= candidate_neighbours.get( nodeAt[nextNode] ).get(candidate).costs;
								}
								saving += NOptList.get(0).oldLength - newLength;
	
								//Is the sum of all previously removed and added edges positive?
								if ( saving > 0.00001 && neighbour_found )
								{
									//remove an existing edge
									int suc = new_neighbour+1;
									int pre  = new_neighbour-1;
									if ( new_neighbour == num_nodes-1)
										suc = 0;
									else if ( new_neighbour == 0)
										pre = num_nodes-1;
									//try to remove the preceding edge of new_neighbour
									if (!containsEdge( NOptList.get(0).brokenEdges, new Edge(new_neighbour, pre)) && pre!=endNode)
									{
										ArrayList<Edge> addedEdges_copy = copyEdgeList(NOptList.get(0).addedEdges);
										ArrayList<Edge> brokenEdges_copy = copyEdgeList(NOptList.get(0).brokenEdges);
										addedEdges_copy.add( new Edge( nextNode, new_neighbour ) );
										brokenEdges_copy.add( new Edge(new_neighbour, pre ) );
										oldLength = model.getDistance( r.nodes.get(new_neighbour).index, r.nodes.get(pre).index );
										degree[new_neighbour]--;
										degree[pre]--;
										NOptMove nextMove = new NOptMove(saving, oldLength, pre, brokenEdges_copy, addedEdges_copy, degree);
										NOptList.add( nextMove );
									}
									//try to remove the succeeding edge of new_neighbour
									if (!containsEdge( NOptList.get(0).brokenEdges, new Edge(new_neighbour, suc)) && suc!=endNode)
									{	
										ArrayList<Edge> addedEdgesC = copyEdgeList(NOptList.get(0).addedEdges);
										ArrayList<Edge> brokenEdgesC = copyEdgeList(NOptList.get(0).brokenEdges);
										addedEdgesC.add( new Edge( nextNode, new_neighbour ) );
										brokenEdgesC.add( new Edge(new_neighbour, suc ) );
										oldLength = model.getDistance(r.nodes.get(new_neighbour).index, r.nodes.get(suc).index);
										degree[new_neighbour]--;
										degree[suc]--;
										NOptMove nextMove = new NOptMove(saving, oldLength, suc, brokenEdgesC, addedEdgesC, degree);
										NOptList.add( nextMove );
									}										
								}//if saving >0
							}//for width
						    NOptList.remove(0);
						}//for list size
					}//move continuation
					if (NOptList.isEmpty())
						{break;}
					else
						depth++;
				}// depth
	
				bestMove = bestNOptMove.copy();		
				//-----------------------------------------------
				//EXECUTE THE BEST MOVE
				//-----------------------------------------------
				if ( bestMove.saving>0.0000001 )
				{
					model.ex1++;
					//System.out.println( "IN  " + bestMove.brokenEdges.size() + "-opt" );
					
					//for each node, get the current neighbours
					ArrayList<ArrayList<Integer>> neighbors = new ArrayList<ArrayList<Integer>>();
					for (int i=0; i<num_nodes; i++)
					{
						ArrayList<Integer> temp = new ArrayList<Integer>();
						if (i==0)
							{temp.add(1); temp.add(num_nodes-1);}
						else if (i==num_nodes-1)
							{temp.add(num_nodes-2); temp.add(0);}
						else
							{temp.add(i-1); temp.add(i+1);}
						neighbors.add(temp);
					}
					//remove neighbours of removed edges
					for (int i=0; i<bestMove.brokenEdges.size(); i++)
					{
						neighbors.get(bestMove.brokenEdges.get(i).node1).remove(new Integer(bestMove.brokenEdges.get(i).node2));
						neighbors.get(bestMove.brokenEdges.get(i).node2).remove(new Integer(bestMove.brokenEdges.get(i).node1));
					}
					//add neighbours of added edges
					for (int i=0; i<bestMove.addedEdges.size(); i++)
					{
						neighbors.get(bestMove.addedEdges.get(i).node1).add(new Integer(bestMove.addedEdges.get(i).node2));
						neighbors.get(bestMove.addedEdges.get(i).node2).add(new Integer(bestMove.addedEdges.get(i).node1));
					}
					//construct the route on the basis of the updated neigbour list
					ArrayList<Node> nodes = new ArrayList<Node>();
					//start with the depot
					nodes.add(r.nodes.get(0));
					nextNode = neighbors.get(0).get(0);
					int prev = 0;
					while (nextNode!=0)//while we have no reached the depot again
					{
						//add one neighbour N
						nodes.add(r.nodes.get(nextNode));
						//determine the other neighbour of N
						if (neighbors.get(nextNode).get(0)!=prev)
							{prev=nextNode; nextNode = neighbors.get(nextNode).get(0);}
						else
							{prev=nextNode; nextNode = neighbors.get(nextNode).get(1);}
					}
					nodes.add(r.nodes.get(0));
					r.nodes = nodes;
					
					//finish
					if (startEdge >= current_edges.size())
						startIndex = startEdge - current_edges.size();
					else
						startIndex = startEdge;
					improved = true;
					
					//keep track of the position of each node within the route
					position = new HashMap<Integer,Integer>();
					for (int i=0; i<num_nodes; i++)
						position.put(r.nodes.get(i).index, i);
					for (int i=0; i<num_nodes; i++)
						nodeAt[i] = r.nodes.get(i).index;
					
					//restart
					break;
				}//execute move	
			}//for start edge
							
		}//while improved
		return r;		
	}	
	
	private static boolean areNeighbours(int node1, int node2, int nodes)
	{
		if (node1+1==node2 || node2+1==node1)
			return true;
		else if(node1==0 && node2==nodes-1)
			return true;
		else if(node2==0 && node1==nodes-1)
			return true;
		else return false;	
	}
	
	private static boolean hasSubTour(ArrayList<Edge> addedEdges, ArrayList<Edge> brokenEdges, int nodeCount)
	{
		ArrayList<ArrayList<Integer>> neighbors = new ArrayList<ArrayList<Integer>>();
		for (int i=0; i<nodeCount; i++)
		{
			ArrayList<Integer> temp = new ArrayList<Integer>();
			if (i==0)
				{temp.add(1); temp.add(nodeCount-1);}
			else if (i==nodeCount-1)
				{temp.add(nodeCount-2); temp.add(0);}
			else
				{temp.add(i-1); temp.add(i+1);}
			neighbors.add(temp);
		}
		//broken edges
		for (int i=0; i<brokenEdges.size(); i++)
		{
			neighbors.get(brokenEdges.get(i).node1).remove(new Integer(brokenEdges.get(i).node2));
			neighbors.get(brokenEdges.get(i).node2).remove(new Integer(brokenEdges.get(i).node1));
		}
		//added edges
		for (int i=0; i<addedEdges.size(); i++)
		{
			neighbors.get(addedEdges.get(i).node1).add(new Integer(addedEdges.get(i).node2));
			neighbors.get(addedEdges.get(i).node2).add(new Integer(addedEdges.get(i).node1));
		}
		//construct the route
		int nextNode = neighbors.get(0).get(0);
		int prev = 0;
		int connectedNodes=1;
		while (nextNode!=0 && connectedNodes < nodeCount)
		{
			if (neighbors.get(nextNode).get(0)!=prev)
				{prev=nextNode; nextNode = neighbors.get(nextNode).get(0);}
			else
				{prev=nextNode; nextNode = neighbors.get(nextNode).get(1);}
			connectedNodes++;
		}
		if (nextNode!=0)//that should not happen
			return true;
		if (connectedNodes < nodeCount)
			return true;
		else
			return false;
	}
	
	public static boolean containsEdge(ArrayList<Edge> list, Edge e)
	{
		for (int i=0;i<list.size(); i++)
			if (list.get(i).equals(e))
				return true;
		return false;
	}
	
	private static ArrayList<Edge> copyEdgeList(ArrayList<Edge> list)
	{
		ArrayList<Edge> copy = new ArrayList<Edge>();
		for (int i=0;i<list.size();i++)
			copy.add( list.get(i).copy() );
		return copy;
	}
	
}
