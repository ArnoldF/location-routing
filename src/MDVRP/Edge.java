package MDVRP;

public class Edge implements Comparable<Edge>
{
	public int node1;
	public int node2;
	public double costs;
	
	//for VRPmodel
	public Edge(int i, int j)
	{
			node1 = i;
			node2 = j;
	}
	
	//for visualization
	public Edge(double c, int i, int j)
	{
			costs = c;
			node1 = i;
			node2 = j;
	}
	
	public Edge copy()
	{
		//return new Edge(this.costs, this.node1,this.node2);
		return new Edge( this.node1,this.node2);
	}
	
	public int compareTo(Edge o) 
	{
		if(o.costs<this.costs)
			return -1;
		else if(o.costs == this.costs)
			return 0;
		else
			return 1;
	}
	
	public boolean equals(Edge e)
	{
		if(e.node1==this.node1 && e.node2==this.node2)
			return true;
		else if(e.node1==this.node2 && e.node2==this.node1)
			return true;
		else
			return false;
	}
	

	
}//class
