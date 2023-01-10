package MDVRP;

import java.util.ArrayList;
import java.util.HashMap;

public class Node 
{
	public int index;
	public int routeIndex;
	public int routePosition;
	public double x;
	public double y;
	public int demand;
	public int[] demandPerProduct;
	public int depot;
	public float detour;//detour induced by the visit of this customer on its route
	public float detour_p;
	public Node[] neighbours;//the two route neighbours
	public boolean isCustomer;
	
	public HashMap<Integer, Float> distances;
	public HashMap<Integer, Float> distances_p;
	public HashMap<Integer, Float> distances_for_CW;	
	public HashMap<Integer, Integer> penalties;	
	
	public HashMap<Integer,InsertionCost>  insert_next_to;	//pre-processed insertion information for RC
	public ArrayList<Node>   is_near_customer_of;
	public ArrayList<Integer>   near_customers;
	public ArrayList<Integer>   customers_considered_as_neigbours;

	
	public Node copy()//DO I EVER NEED THIS? DO I EVER WANT TO RESET?
	{
		Node n	= new Node(index);
		n.routeIndex	= routeIndex;
		n.routePosition	= routePosition;
		n.x				= x;
		n.y				= y;
		n.demand		= demand;
		n.isCustomer    = isCustomer;
		n.demandPerProduct		= demandPerProduct.clone();
		//n.demandPerProduct = demandPerProduct.clone();
		if(distances!=null) n.distances		= (HashMap<Integer, Float>) distances.clone();
		//if(penalties!=null) n.penalties     = new int[penalties.length];//TODO create a hashmap based on distances with all values set to 0
		if(distances_p!=null) n.distances_p	= (HashMap<Integer, Float>) distances_p.clone();//distances //TODO does this work?
		if(distances_for_CW!=null) n.distances_for_CW	= (HashMap<Integer, Float>) distances_for_CW.clone();
		if(penalties!=null) n.penalties	= (HashMap<Integer, Integer>) penalties.clone();
		if(near_customers!=null) n.near_customers	= (ArrayList<Integer>) near_customers.clone();
		if(customers_considered_as_neigbours!=null) n.customers_considered_as_neigbours	= (ArrayList<Integer>) customers_considered_as_neigbours.clone();
		//if(isNeighbourOf!=null) n.isNeighbourOf	= (ArrayList<Node>) isNeighbourOf.clone();
		n.depot			= depot;
		n.detour_p		= detour_p;
		if(neighbours!=null) n.neighbours    = neighbours.clone();
		//TODO create class copes in a map
		//if(nearestNeighbours!=null) 
			//n.nearestNeighbours =  new HashMap<Integer, InsertionCost>(nearestNeighbours);
		//if(nearestNeighboursT!=null) 
			//n.nearestNeighboursT = (HashMap<Integer, InsertionCost>) nearestNeighboursT.clone();
		/*for (int i=0; i<nearestNeighbours..size(); i++ )
			//if (nearestNeighboursT[i]!=null)
			n.nearestNeighbours.put(nearestNeighbours., value) = nearestNeighbours[i].copy();
		if(nearestNeighboursT!=null) 
			for (int i=0; i<nearestNeighboursT.length; i++ )
				if (nearestNeighboursT[i]!=null)
					n.nearestNeighboursT[i] = nearestNeighboursT[i].copy();*/
				
		n.detour		= detour;

		return n;
	}
	
	
	public ArrayList<Node> getNeighbours() {
		ArrayList<Node> neighbours = new ArrayList<>();
		if (this.neighbours[0].isCustomer)
			neighbours.add(this.neighbours[0]);
		if (this.neighbours[1].isCustomer)
			neighbours.add(this.neighbours[1]);
		return neighbours;
		
	}
	
	public Node(int i)
	{
			index = i;
			distances			= new HashMap<Integer, Float>();
			distances_p			= new HashMap<Integer, Float>();
			distances_for_CW	= new HashMap<Integer, Float>();
			penalties			= new HashMap<Integer, Integer>();
			insert_next_to 	= new HashMap<Integer, InsertionCost>();
			is_near_customer_of      	= new ArrayList<Node>();
			near_customers = new ArrayList<Integer>();
			customers_considered_as_neigbours = new ArrayList<Integer>();
			neighbours = new Node[2];
	}
	
}//class
