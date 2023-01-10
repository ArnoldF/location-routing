package MDVRP;

import java.util.ArrayList;

public class Route 
{
	public ArrayList<Node> nodes;
	public int length; //number of customers on the route
	public int load;
	public double cost;
	public int depot;
	
	public Route()
	{
		 nodes = new ArrayList<Node>();
		 this.cost	= 0;
		 this.load	= 0;
		 this.length	= 0;
	}
		
	public void add(Node n)
	{
		nodes.add(n);
	}
	
	public void clear()
	{
		this.cost	= 0;
		this.load	= 0;
		this.length	= 0;
		this.depot  = -1;
		this.nodes.clear();	
	}
	
	public Route copy()
	{
		Route r 	= new Route();
		r.cost 		= this.cost;
		r.load 		= this.load;
		r.length 	= this.length;
		r.depot	    = this.depot;
		r.nodes 	= new ArrayList<Node>();
		for (int i=0; i<this.nodes.size(); i++)
			r.nodes.add(this.nodes.get(i));
		return r;
	}
	
	public Route copy(Node[] nodesNew)
	{
		Route r 	= new Route();
		r.cost 		= this.cost;
		r.load 		= this.load;
		r.length 	= this.length;
		r.nodes 	= new ArrayList<Node>();
		r.depot     = this.depot;
		for (int i=0; i<this.nodes.size(); i++)
			r.nodes.add( nodesNew[this.nodes.get(i).index] );
		return r;
	}
	
	public void print()
	{
		System.out.print(cost + ": ");
		for (int i=0; i<length+2; i++)
			System.out.print(nodes.get(i) + ", ");
		System.out.println();
	}

}//class
