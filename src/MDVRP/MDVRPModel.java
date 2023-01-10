package MDVRP;

import java.io.File;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;



public class MDVRPModel {

	//model parameters
	public  int customers;	
	public int depots;
	public int vertices;
	public int products;
	public int[][] depotInventory;
	public int vehicles;
	public int capacityLimitVehicle;
	public double maxTourLength;
	public boolean distances_rounded;
	public double cost_per_Route  = 0;
	
	//run parameters			
	public double baselineUnit;
	private double wLength = 1.0;
	private double wWidth = 1.0;
	public boolean penalties_off;
	public boolean inventoryOff = true;
	public boolean heuristicStart = true;
	double penalty = 0;
	int perturbation_length;
	
	//model variables
	public Node[] nodes;
	public ArrayList<Route> routes;
	public ArrayList<Integer> routes_to_be_optimised;
	public ArrayList<Integer> routes_update_edge_values;
	public ArrayList<Integer> routes_changed_during_perturbation;
	public ArrayList<Integer> routes_update_RC_variables;
	public ArrayList<Integer> starts_RC;
	public ArrayList<Edge>	starts_CROSS;
	
	private Edge[][] edgeValues;
	public int [][] remainingInventory;
	public int moves_RC = 0;
	public int moves_CROSS = 0;
	int penalties_given = 0;
	ArrayList<EdgeClassification> edgeClasses;
	public Edge[] distance_to_nearest_depot;
	
	double[] development;
	int developmentCount = 0;
	int ev1 = 0;
	int ev2 = 0;
	int ev3 = 0;
	int ex1 = 0;
	int ex2 = 0;
	int ex3 = 0;
	int minNumberRoutes;
	
	//MTVRP
	double maxTourLength2 = 99999;
	double dayLimit = 0;
	int maxVehicles = 0;
	
	//LRP
	Node[] possibleDepots;
	Node hub;
	int capacityLimitVehicleFirstLevel;
	boolean isFirstLevel = false;
	
	public int pruneLimit = 30;
	
	//pattern Mining
	public HashMap<Long, Integer> minedPatterns;
	ArrayList<ArrayList<Long>> patternHistory;
	public HashMap<Long, Integer> mostFrequPattern = new HashMap<Long, Integer>();
	int numConsideredPattern = 1000;
	int minMaxFrequency = 0;
	
	int patternCollected = 0;
	int statsCollected = 0;
	int[] averageStats = new int[11];
	int[] averageStatsMoves = new int[11];
	int[][] averageMoveChanges = new int[10][11];
	int[][] averageRouteChanges = new int[10][11];
	int[] averageNumNodes = new int[11];
	int maxCollection = 10;
	int maxPatternSize = 6;//TODO if this is too large then the hash key becomes to large
	int minPatternSize = 3;
	int countPatternImprove;
	int countPatternTotal;
	int countPatternNotCurrent;
	
	LinKernighan lk;
	
	public static class EdgeClassification
	{		
		double wWidth;
		double wLength;
		
		public EdgeClassification(double width, double length)
		{
			wWidth		= width;
			wLength		= length;
		}
		public EdgeClassification copy()
		{
			return new EdgeClassification(this.wWidth, this.wLength);
		}
	}
	
	public static class EdgeSequence
	{		
		int accCapacity;
		int startNode;
		int endNode;
		ArrayList<Integer> nodes;
		
		public EdgeSequence(int a, int s, int e, ArrayList<Integer> n)
		{
			accCapacity = a;
			startNode = s;
			endNode = e;
			nodes = n;
		}
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//FUNCTION : Constructor
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	public MDVRPModel(int customers, int depots, int products)
	{
		this.customers 		= customers;
		this.depots    		= depots;
		this.vertices  		= customers + depots;
		nodes 				= new Node[vertices];
		maxTourLength 		= Float.MAX_VALUE;
		//inventory
		remainingInventory 	= new int[depots][products];
		depotInventory 	   	= new int[depots][products];
		this.products  		= products;
		//lists
		routes 						= new ArrayList<Route>();
		routes_to_be_optimised 		= new ArrayList<Integer>();
		routes_changed_during_perturbation = new ArrayList<Integer>();
		starts_RC 					= new ArrayList<Integer>();
		routes_update_RC_variables 	= new ArrayList<Integer>();
		routes_update_edge_values 	= new ArrayList<Integer>();
		starts_CROSS     			= new ArrayList<Edge>();
		distance_to_nearest_depot  	=  new Edge[customers];
	}
	
	
	public MDVRPModel copy() {
		MDVRPModel copiedModel = new MDVRPModel(this.customers, this.depots, this.products);
		
		copiedModel.capacityLimitVehicle = this.capacityLimitVehicle;
		copiedModel.cost_per_Route = this.cost_per_Route;
		copiedModel.nodes 	= new Node[customers];
		copiedModel.perturbation_length = 30;
		copiedModel.inventoryOff = false;
		copiedModel.distances_rounded = false;
		copiedModel.pruneLimit = Math.min(30, customers-1);
		
		for (int i=0; i<customers; i++) {
			copiedModel.nodes[i] = this.nodes[i].copy();			
		}
		for (int i=0; i<customers; i++) {
			for (int j=0; j<pruneLimit; j++) {
				copiedModel.nodes[i].insert_next_to.put(copiedModel.nodes[i].near_customers.get(j), new InsertionCost(0,0,-1,-1));
				copiedModel.nodes[copiedModel.nodes[i].near_customers.get(j)].is_near_customer_of.add(copiedModel.nodes[i]);
			}
		}
		

		//lists
		copiedModel.routes 						= new ArrayList<Route>();
		copiedModel.routes_to_be_optimised 		= new ArrayList<Integer>();
		copiedModel.routes_changed_during_perturbation = new ArrayList<Integer>();	
		copiedModel.starts_RC 					= new ArrayList<Integer>();
		copiedModel.routes_update_RC_variables 	= new ArrayList<Integer>();
		copiedModel.routes_update_edge_values 	= new ArrayList<Integer>();
		copiedModel.starts_CROSS     			= new ArrayList<Edge>();
		copiedModel.distance_to_nearest_depot  	=  new Edge[customers];
		copiedModel.maxTourLength = Double.MAX_VALUE;
		
		return copiedModel;
	}
	
	
	
	public float getDistance(int n1, int n2)
	{
		if (nodes[n1].distances.containsKey(n2))
			return nodes[n1].distances.get(n2);
		else if (nodes[n2].distances.containsKey(n1))
			return nodes[n2].distances.get(n1);
		else
		{
			//System.out.println("A requested distance value is not stored");
			float dist = 0;
			if (distances_rounded)
				dist = 100 * (int) (0.5 + Math.sqrt( Math.pow( nodes[n1].x - nodes[n2].x , 2) + Math.pow( nodes[n1].y - nodes[n2].y , 2) ) );
			else
				dist = 100 * (float) Math.sqrt( Math.pow( nodes[n1].x - nodes[n2].x , 2) + Math.pow( nodes[n1].y - nodes[n2].y , 2) );
			return dist;
		}
	}
	
	public float getDistance2(int n1, int n2)
	{
		if (n1<customers && nodes[n1].distances_p.containsKey(n2))
			return nodes[n1].distances_p.get(n2);
		else if (n2<customers && nodes[n2].distances_p.containsKey(n1))
			return nodes[n2].distances_p.get(n1);
		else
		{
			//System.out.println("A requested distance value is not stored");
			int n1t = Math.min(n1, customers);
			int n2t = Math.min(n2, customers);
			float dist = 0;
			if (distances_rounded)
				dist = (int) (0.5 + Math.sqrt( Math.pow( nodes[n1t].x - nodes[n2t].x , 2) + Math.pow( nodes[n1t].y - nodes[n2t].y , 2) ));
			else
				dist = (float) Math.sqrt( Math.pow( nodes[n1t].x - nodes[n2t].x , 2) + Math.pow( nodes[n1t].y - nodes[n2t].y , 2) );			
			dist += penalty * getPenalties(n1, n2);
			return dist;
		}
	}
	
	public void setDistance2(int n1, int n2, float val)
	{
		if (n1< customers && nodes[n1].distances_p.containsKey(n2))
			nodes[n1].distances_p.put(n2, val);
		if (n2< customers && nodes[n2].distances_p.containsKey(n1))
			nodes[n2].distances_p.put(n1, val);
		/*else
		{
			System.out.println("A distance value cannot be set");
		}*/
	}
	
	public int getPenalties(int n1, int n2)
	{
		if (n1< customers && nodes[n1].penalties.containsKey(n2))
			return nodes[n1].penalties.get(n2);
		else if (n2< customers && nodes[n2].penalties.containsKey(n1))
			return nodes[n2].penalties.get(n1);
		else
		{
			//System.out.println("A requested distance value is not stored");
			return 0;
		}
	}
	
	
	public void setPenalties(int n1, int n2, int val)
	{
		if (n1< customers && nodes[n1].penalties.containsKey(n2))
			nodes[n1].penalties.put(n2, val);
		else
			if (n1 < customers)
				nodes[n1].penalties.put(n2, val);
		if (n2< customers && nodes[n2].penalties.containsKey(n1))
			nodes[n2].penalties.put(n1, val);
		else
			if (n2 < customers)
				nodes[n2].penalties.put(n1, val);
	}
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//FUNCTION : Compute the costs of the model
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	public double computeSolutionCosts() {
		double totalCost = 0;
		for (int r=0; r<routes.size(); r++)
			if (routes.get(r).length > 0)
			for (int i=1; i<routes.get(r).length+2; i++)
				totalCost += getDistance(routes.get(r).nodes.get(i-1).index, routes.get(r).nodes.get(i).index);
		
		//LRP
		for (int i=0; i<routes.size(); i++ )
			if ( routes.get(i).length > 0)
				totalCost+=cost_per_Route;
		return totalCost;
	}
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//FUNCTION : Check feasibility of solution
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	public boolean checkSolution()
	{
		boolean validSolution = true;
		int[][] inventoryUsed = new int[depots][Math.max(1,products)];
		boolean visited[] = new boolean[customers+1];
		for (int r=0; r<routes.size(); r++)
		{
			int load = 0;
			for (int i=1; i<routes.get(r).nodes.size()-1; i++)
			{
				load 		+= routes.get(r).nodes.get(i).demand;
				if (visited[routes.get(r).nodes.get(i).index])
					System.out.println("ERROR : Customer assigned twice - " + routes.get(r).nodes.get(i).index);
				visited[routes.get(r).nodes.get(i).index]	=	true;
				if ( routes.get(r).nodes.get(i).routePosition != i  )
					System.out.println("ERROR : route position wrong " + routes.get(r).nodes.get(i).index);
				if ( routes.get(r).nodes.get(i).routeIndex != r  )
					System.out.println("ERROR : route index wrong " + routes.get(r).nodes.get(i).index);
				for (int p=0; p<Math.max(1,products); p++)
					inventoryUsed[routes.get(r).depot][p] += routes.get(r).nodes.get(i).demandPerProduct[p];
			}
			//check cost computation and tour length constraints
			double cost = 0;
			for	(int i=1; i<routes.get(r).length+2; i++)
				 	cost += getDistance( routes.get(r).nodes.get(i-1).index, routes.get(r).nodes.get(i).index);
			if ( Math.abs(cost- routes.get(r).cost) > 0.01)
				;//System.out.println("ERROR : wrongCost "+r+";"+" "+routes.get(r).cost +" instead of "+ cost);
			if ( cost > this.maxTourLength)
			{
				validSolution = false;
				//System.out.println("ERROR : cost too high on route " + r);
			}
			//check capacity constraints
			if (load > capacityLimitVehicle && routes.get(r).length>1)//large
			{
				validSolution = false;
				System.out.println("ERROR : capacity violation on route " + r);
			}
			if (load != routes.get(r).load)
				System.out.println("ERROR : load computation wrong");
			if (routes.get(r).nodes.get( routes.get(r).nodes.size()-1 ).index < customers)
				System.out.println("ERROR : last node is not a depot - " + r);		
			if (routes.get(r).nodes.get( 0 ).index != routes.get(r).nodes.get( routes.get(r).nodes.size()-1 ).index)
				System.out.println("ERROR : route starts and ends at different depots");
			if (routes.get(r).nodes.get( 0 ).index-customers != routes.get(r).depot)
				System.out.println("ERROR : route has wrong depot index");
		}
		//check whether all customers have been visited
		for (int i=0; i<customers; i++)
			if (!visited[i] && nodes[i]!=null)
				System.out.println("ERROR : customer " + i + " not visited");
		
		//check inventory constraints
		//if (!this.inventoryOff)
		//TODO
		for (int d=0; d<depots; d++)
			for (int p=0; p<Math.max(1,products); p++)
			{
				if (inventoryUsed[d][p] > this.depotInventory[d][p])
				{
					;//System.out.println("ERROR : inventory overused");
					validSolution = false;
				}
				//if (this.remainingInventory[d][p] + inventoryUsed[d][p] != this.depotInventory[d][p])
					//System.out.println("inventory comp wrong");
			}
				
		return validSolution;
	}
	
	public boolean checkINventoryConstraints() {
		int[] inventoryUsed = new int[depots];
		for (int r=0; r<routes.size(); r++)
			for (int i=1; i<routes.get(r).nodes.size()-1; i++)
				inventoryUsed[routes.get(r).depot] += routes.get(r).nodes.get(i).demand;
		
		for (int d=0; d<depots; d++) {
			if (inventoryUsed[d] > this.depotInventory[d][0]) {
				//System.out.println("ERROR : inventory overused");
				return false;
			}
		}
		return true;
	}
	
	
	public void intraRouteOptimization(int max, boolean optimise)
	{
		for ( int i=0; i<routes.size(); i++ )
		{
			if ( routes_to_be_optimised.contains(i) )
			{
				if (optimise)
				{
					Route r = lk.linKernighan(routes.get(i), max, this);
					//replace old route with optimised route
					routes.remove(i);
					routes.add(i, r);
					for (int j=1; j<routes.get(i).length+1; j++)
					{
						routes.get(i).nodes.get(j).routePosition = j;
						routes.get(i).nodes.get(j).routeIndex = i;
					}				
				}
				//compute the new costs
				double cost = 0;
				for (int j=1; j<routes.get(i).length+2; j++)
				 	cost += this.getDistance(routes.get(i).nodes.get(j-1).index, routes.get(i).nodes.get(j).index);//routes.get(i).nodes.get(j-1).distances[routes.get(i).nodes.get(j).index];
				routes.get(i).cost = cost;
				//update the values for each edge later
				if (!routes_update_edge_values.contains(i)) 
					routes_update_edge_values.add( i );
			}
			
		}
		update_neighbours();
		if (optimise)
			routes_to_be_optimised.clear();
	}
	
	//after changing the order of a route, update the neighbours of each customer
	//note that we have to do this for the whole route, since direction of subtours might change
	public void update_neighbours()
	{
		for (int r=0; r<routes_to_be_optimised.size(); r++)
		{
			int route = routes_to_be_optimised.get(r);
			for (int n=1; n<routes.get(route).length+1; n++)
			{
				Node node = routes.get(route).nodes.get(n);
				node.neighbours[0] = routes.get(node.routeIndex).nodes.get(node.routePosition - 1);
				node.neighbours[1] = routes.get(node.routeIndex).nodes.get(node.routePosition + 1);
			}
		}
	}
		
	public boolean hasSufficientInventory(int depot, int[] demandPerProduct)
	{
		boolean sufficient = true;
		for (int p=0; p<products; p++)
			if ( this.remainingInventory[depot][p] < demandPerProduct[p])
			{
				sufficient = false;
				break;
			}
		return sufficient;
	}
		
	public boolean haveSufficientInventory(int depot1, int depot2, int[] demandPerProduct1, int[] demandPerProduct2)
	{
		boolean sufficient = true;
		for (int p=0; p<products; p++)
			if ( this.remainingInventory[depot1][p] < demandPerProduct2[p] - demandPerProduct1[p] ||
			     this.remainingInventory[depot2][p] < demandPerProduct1[p] - demandPerProduct2[p])
			{
				sufficient = false;
				break;
			}
		return sufficient;
	}
		
	//*******************************************************************************
	//The following two functions update the pre-processed data for the Relocation chain
	//1) the detour costs of each node
	//2) the insertion costs next to a nearby node
	//*******************************************************************************
	public void update_parameters_for_EC()
	{
		for ( int r=0; r < routes_update_RC_variables.size(); r++ )
		{
			int route = routes_update_RC_variables.get( r );
			for (int n=1; n<routes.get( route ).length+1; n++)
			{
				Node node = routes.get( route ).nodes.get(n);
				if ( node.index < customers )
				{
					//1.) update the detour costs in terms of penalized edge costs (detour) and the true detour costs (detourT)
					node.detour_p =  getDistance2(node.index, node.neighbours[0].index);
					node.detour_p += getDistance2(node.index, node.neighbours[1].index);
					node.detour_p -= getDistance(node.neighbours[0].index, node.neighbours[1].index);
					
					node.detour = getDistance(node.index, node.neighbours[0].index);
					node.detour+= getDistance(node.index, node.neighbours[1].index);
					node.detour-= getDistance(node.neighbours[0].index, node.neighbours[1].index);
					
					//2.) update the insertion costs next to nearest neighbours
					for (int nearest_node=0; nearest_node<node.is_near_customer_of.size(); nearest_node++)
						if (nodes[node.is_near_customer_of.get(nearest_node).index] != null)
						{
							updateNearestNeighbour(node.is_near_customer_of.get(nearest_node), node);		
						}
				}
			}
		}
		routes_update_RC_variables.clear();
	}		
				
	private void updateNearestNeighbour(Node node, Node neighbour)
	{
		if ( node.routeIndex != neighbour.routeIndex )
		{
			int n1 = neighbour.neighbours[0].index;
			int n2 = neighbour.neighbours[1].index;
			double dist1;
			double dist2;
			double dist1T;
			double dist2T;

			dist1 = getDistance(node.index, neighbour.index) + getDistance(node.index, n1) - getDistance2(n1, neighbour.index);//node.distances[neighbour.index] + node.distances[n1] - neighbour.distances2[n1];
			dist2 = getDistance(node.index, neighbour.index) + getDistance(node.index, n2) - getDistance2(n2, neighbour.index);//node.distances[neighbour.index] + node.distances[n2] - neighbour.distances2[n2];
			dist1T = getDistance(node.index, neighbour.index) + getDistance(node.index, n1) - getDistance(n1, neighbour.index);//node.distances[neighbour.index] + node.distances[n1] - neighbour.distances[n1];
			dist2T = getDistance(node.index, neighbour.index) + getDistance(node.index, n2) - getDistance(n2, neighbour.index);//node.distances[neighbour.index] + node.distances[n2] - neighbour.distances[n2];
					
			if (dist1 < dist2)
				node.insert_next_to.put(neighbour.index, new InsertionCost(dist1, dist1T, neighbour.routeIndex, neighbour.routePosition));
			else
				node.insert_next_to.put(neighbour.index, new InsertionCost(dist2, dist2T, neighbour.routeIndex, neighbour.routePosition+1));

		}
		else //nodes are already in the same route
			node.insert_next_to.put(neighbour.index, new InsertionCost(Double.MAX_VALUE, Double.MAX_VALUE, -1, -1));
	}
	
	
	public double computeEdgeWidth(double centerX, double centerY, Node node1, Node node2, int route)
	{
		//compute the distance to the center, for both nodes
		Node depot 	  = routes.get(route).nodes.get(0);
		double dist1  = (centerY - depot.y) * node1.x;
		dist1 		 -= (centerX - depot.x) * node1.y;
		dist1		 += (centerX * depot.y) - (centerY * depot.x);//TODO this is always the same for a route
		dist1         = dist1 / Math.sqrt( Math.pow(centerY - depot.y, 2) + Math.pow(centerX - depot.x, 2) );
		
		double dist2  = (centerY - depot.y) * node2.x;
		dist2 		 -= (centerX - depot.x) * node2.y;
		dist2		 += (centerX * depot.y) - (centerY * depot.x);
		dist2         = dist2 / Math.sqrt( Math.pow(centerY - depot.y, 2) + Math.pow(centerX - depot.x, 2) );
		
		double width = Math.abs( dist1 - dist2 );
		return width;
	}
	
	//*******************************************************************************
	//The following functions initialize and update the "badness" of an edge
	//*******************************************************************************
	private void initEdgeValues()
	{
		edgeValues = new Edge[customers][2];
		for ( int n = 0; n < customers; n++)
		{
			edgeValues[n][0]  	= new Edge(0, -1, -1);
			edgeValues[n][1]  	= new Edge(0, -1, -1);
		}
	}
	
	private void updateEdgeValues()
	{
		for (int i=0; i < routes_update_edge_values.size(); i++)
		{
			int route 	= routes_update_edge_values.get( i );
			double centerX = computeRouteCenterX(route);
			double centerY = computeRouteCenterY(route);
			
			for (int j=1; j<routes.get(route).length+2; j++)
			{			
				int node1 	= routes.get(route).nodes.get(j-1).index;
				int node2   = routes.get(route).nodes.get(j).index;

				//compute the "badness" of this edge
				double cost = 0;
				//length
				cost += wLength * getDistance(node1, node2);
				//width
				cost += wWidth * computeEdgeWidth(centerX, centerY, nodes[node1], nodes[node2], route);
				//divide by previously received penalties
				if (nodes[node1].penalties.containsKey(node2) || nodes[node2].penalties.containsKey(node1))
					cost = cost / (1.0 + this.getPenalties(node1, node2));
				
				//assign edge value to both customers
				if (node1 < customers)
					edgeValues[node1][1]  	= new Edge(cost, node1, node2);
				if (node2 < customers)
					edgeValues[node2][0]  	= new Edge(cost, node2, node1);
			}
		}
		routes_update_edge_values.clear();
	}
	
	public double computeRouteCenterX(int route)
	{
		double xCenter = 0;
		for (int j=0; j< routes.get(route).nodes.size()-1; j++)
			xCenter += routes.get(route).nodes.get(j).x;
		return xCenter / (double)(routes.get(route).nodes.size()-1);
	}
	
	public double computeRouteCenterY(int route)
	{
		double yCenter = 0;
		for (int j=0; j< routes.get(route).nodes.size()-1; j++)
			yCenter += routes.get(route).nodes.get(j).y;
		return yCenter / (double)(routes.get(route).nodes.size()-1);
	}
	
	
	public void optimizeRoutes(int maxIterations, double maxRunTime, boolean messagesOn)//boolean compare;  double limit - both for data mining study
	{		
		//******************************************************************
		//VARIABLE INITIALISATION
		//******************************************************************
		//penalization criterion
		edgeClasses = new ArrayList<EdgeClassification>();
		//width
		edgeClasses.add( new EdgeClassification(1,0) );
		//length
		edgeClasses.add( new EdgeClassification(0,1) );
		//width plus length
		edgeClasses.add( new EdgeClassification(1,1) );
		
		//parameters
		int lambda 				= 4;	
		
		//stopping and measuring time
		Stopwatch stopwatch = new Stopwatch();
		stopwatch.start();
		Stopwatch stopwatchOperators = new Stopwatch();
		double run_time_operators[] = new double[4];
		double run_time_to_Best = 0;
		double run_time_to_Restart = 0;
		double run_time = 0;
		
		EjectionChain ejectionChain = new EjectionChain();
		CROSS cross = new CROSS();
		lk = new LinKernighan();
		
		minedPatterns = new HashMap<Long, Integer>();
		patternHistory = new ArrayList<ArrayList<Long>>();
		
		//save progress
		developmentCount = 1;
		development = new double[1500];
		
		//save global optimum
		ArrayList<Route> best_routing = new ArrayList<Route>();
		int best_iteration = 0;
		
		//variables to control the flow of the heuristic	
		heuristicStart 				= true;
		int iteration_count 		= 0;
	    moves_RC 					= 0;
	    moves_CROSS 				= 0;
	    ex1 = ex2 = ex3 =0;
	    int edgeClassIndex			= 0;
	    ArrayList<HashMap<Integer, Float>> distances2_copy = new ArrayList<HashMap<Integer, Float>>();
		
	    //edge penalisation
	    wWidth	= 1;
	    wLength	= 0;
	    this.penalties_off = false;

		
		//initialize penalized distances
		for (int i=0; i<nodes.length; i++)
			nodes[i].distances_p   = (HashMap<Integer, Float>) nodes[i].distances.clone();
		
		//initialize lists
		routes_update_RC_variables.clear();
		routes_update_edge_values.clear();
		routes_to_be_optimised.clear();
		routes_changed_during_perturbation.clear();
		starts_CROSS.clear();
		ArrayList<Edge> penalizedEdges = new ArrayList<Edge>();
		for(int i=0; i< routes.size(); i++)
		{
			routes_to_be_optimised.add(i);
			routes_update_edge_values.add(i);
			routes_update_RC_variables.add( i );
		}
		
	    //only MDVRP: distance to closest depot
		    for(int i=0; i<customers; i++)
			{
				double minDist = 999999;
				int nearestDepot = -1;
				for(int j=0; j<depots; j++)
					if ( this.getDistance(i, customers+j) < minDist )
					{
						minDist = this.getDistance(i, customers+j);
						nearestDepot = j;
					}
				distance_to_nearest_depot[i] = new Edge( minDist, nearestDepot, -1 );
			}
		    
		//******************************************************************
		//HEURISTIC INITIALISATION
		//******************************************************************	
		double best_solution_value = computeSolutionCosts();
		double best_solution_MTVRP = 9999999;
	    development[0] = best_solution_value;
	    //update applet
		
		
		//intra-route optimisation
		intraRouteOptimization(lambda, true);
		
		stopwatch.stop();
		 run_time += stopwatch.elapsedTimeMillis();
		 stopwatch.start();
		 if (run_time > 999000 * developmentCount)
		 {
			 development[developmentCount] = best_solution_value;
			 developmentCount++;
		 }

		baselineUnit = (float)computeSolutionCosts() / ((float)customers * 10.0);
		penalty 		= baselineUnit;
		initEdgeValues();
		updateEdgeValues();

	 	//safe initial solution
		this.checkSolution();
		for (int r=0; r<routes.size(); r++)
			best_routing.add( routes.get(r).copy() );// copy routes, with references to nodes (adapt routeIndex and routePosition)
	    best_solution_value = computeSolutionCosts();	    
	    
	    //start with an optimisation of the entire solution, RC starts with all customers and CROSS with all edges
	    starts_RC.clear();
		for (int i=0; i<customers; i++)
		    starts_RC.add(i);
		for (int i=0; i<routes.size(); i++)
			for (int j=1; j<routes.get(i).length; j++)
				starts_CROSS.add( new Edge(routes.get(i).nodes.get(j-1).index, routes.get(i).nodes.get(j).index));

		//******************************************************************
		//MAIN LOOP
		//******************************************************************
		outerLoop:
	    while ( iteration_count < maxIterations && run_time < maxRunTime)
	    {
	    	//******************************************************************
			//LOCAL SEARCH
	    	//******************************************************************
			 double previousSolution	= 0;
			 double currentSolution  	= 0;
			 boolean improved1 			= true;
			 
			 stopwatchOperators.start();
			 if(penalties_off)
				 intraRouteOptimization(lambda, penalties_off);
			 stopwatchOperators.stop();
			 run_time_operators[2] += stopwatchOperators.elapsedTimeMillis();
			
			 while (improved1)
			 {
				previousSolution	  	= computeSolutionCosts();	
					 
				//CROSSOVER
				stopwatchOperators.start();		
				cross.cross(this, starts_CROSS);
				stopwatchOperators.stop();
				run_time_operators[1] += stopwatchOperators.elapsedTimeMillis();

				//LIN-KERNIGHAN
				stopwatchOperators.start();
				intraRouteOptimization(lambda,this.penalties_off || heuristicStart);
				stopwatchOperators.stop();
				run_time_operators[2] += stopwatchOperators.elapsedTimeMillis();
					 
				//RELOCATION CHAIN
				stopwatchOperators.start();
				update_parameters_for_EC();		
				//if ( heuristicStart)
					//sort_insertions();
				ejectionChain.ejectionChain(this, 3, starts_RC);					 
				stopwatchOperators.stop();
				run_time_operators[0] += stopwatchOperators.elapsedTimeMillis();
					
				//LIN-KERNIGHAN
				stopwatchOperators.start();
				intraRouteOptimization(lambda,this.penalties_off || heuristicStart);	
				stopwatchOperators.stop();
				run_time_operators[2] += stopwatchOperators.elapsedTimeMillis();
				 		 
				 //EVALUATION
				 currentSolution = computeSolutionCosts();
				 if ( currentSolution < previousSolution && (this.penalties_off || heuristicStart) )
					 previousSolution = currentSolution;
				 else
					 improved1 = false;
				 
				 stopwatch.stop();
				 run_time += stopwatch.elapsedTimeMillis();
				 stopwatch.start();
				 if (run_time > 999000 * developmentCount)
				 {
					 development[developmentCount] = Math.min(currentSolution, best_solution_value);
					 developmentCount++;
				 }
				 
				 if (run_time > maxRunTime)
					 break;
				 
				 //if (iteration_count == 0 && this.moves_CROSS + this.moves_RC == 0)
					 //break outerLoop;
			}
			//---------------------------------------------------------
			//finalize local search
			//---------------------------------------------------------
			if  (this.penalties_off || heuristicStart)
			{
				iteration_count++;
				
				stopwatch.stop();
				 stopwatch.start();
				 
				 
				/*if (patternCollected == maxCollection)
					{extractPattern(minPatternSize, maxPatternSize, true);}
				else{
					extractPattern(minPatternSize, maxPatternSize, false);
					patternCollected++;
				}
				currentSolution = computeSolutionCosts();*/
				
				stopwatch.stop();
				run_time_operators[3] += stopwatch.elapsedTimeMillis();
				 stopwatch.start();
				 
				//if (patternCollected == maxCollection)
					//getPatternStats(maxCollection);
			}
				//sort_insertions();

			 heuristicStart = false;						 
			 starts_RC.clear();
			 starts_CROSS.clear();
			 this.checkSolution();
							
			 //NEW GLOBAL OPTIMUM?
			 if (best_solution_value > currentSolution * 1.00001 && this.checkSolution())
			 {
				 //safe this solution
				 best_solution_value 		= currentSolution;
				 best_iteration 			= iteration_count;
				 best_routing.clear();
				 for (int r=0; r<routes.size(); r++)
					best_routing.add( routes.get(r).copy() );		
				 
				 run_time_to_Best = run_time;		
				 run_time_to_Restart = run_time;
				 
				 int countActiveRoutes = 0;
				 for (int r=0; r<routes.size(); r++)
					 if (routes.get(r).length > 0)
						 countActiveRoutes++;
				if (messagesOn)
				{
					int totalMoves = ex1+ex2+ex3;
					int intraMoves = ex2+ex3;
					System.out.println("Time : " + String.format("%.1f",run_time_to_Best/1000) + "  \t" + "Value : " + String.format("%.2f",best_solution_value) + " (" + countActiveRoutes + " Routes) \t" + "Moves : " + totalMoves );
					//System.out.println("computation times : " + run_time_operators[0] + ", " + run_time_operators[1] + ", " +run_time_operators[2] + ", " +run_time_operators[3] + ": " + countActiveRoutes);
				}
				 
			 }//best solution
			 
			
			 
			 //update run statistics
			 stopwatch.stop();
			 run_time += stopwatch.elapsedTimeMillis();
			 stopwatch.start();		 
			 
			/*if (isFeasibleMTVRP(this.maxVehicles, this.dayLimit))
				if (best_solution_MTVRP > currentSolution)
				{
					System.out.println("OK " + currentSolution + " " + run_time);
					this.checkSolution();
					best_solution_MTVRP = currentSolution;
				}*/


			//------------------------------------------------------------------------
			//CHANGE FROM OPTIMISATION TO PERTURBATION
			//------------------------------------------------------------------------
			//reset penalized edge values after reoptimisation
			if (penalties_off)
			{
				penalties_off = false;
				routes_changed_during_perturbation.clear();

				//copy the penalised distance values
				for (int i=0; i<customers; i++)
					nodes[i].distances_p = (HashMap<Integer, Float>) distances2_copy.get(i).clone();
				for(int i=0; i< routes.size(); i++)
					routes_update_RC_variables.add(i);
				
				edgeClassIndex++;
				 if (edgeClassIndex >= edgeClasses.size())
					 edgeClassIndex = 0;
				 wWidth		= edgeClasses.get(edgeClassIndex).wWidth;
				 wLength	= edgeClasses.get(edgeClassIndex).wLength;

				 //update according to new penalization criterion
				 for(int i=0; i< routes.size(); i++)
					routes_update_edge_values.add(i);
				 
				 /*if ( run_time - run_time_to_Restart > customers * 600 )
				 {
					 resetSolution( best_routing, null, null);
					 for(int i=0; i< routes.size(); i++)
							routes_with_new_edges.add(i);
					 //System.out.println("    RESTART " + wWidth + " "+ wLength + " "+ wDepth + " " + pruneLimit);
					 run_time_to_Restart = run_time;
				 }*/			
			}
			
			//------------------------------------------------------------------------
			//CHANGE FROM PERTURBATION TO OPTIMISATION
			//------------------------------------------------------------------------
			
			if(  ( moves_RC + moves_CROSS >= iteration_count * perturbation_length ) )
			{	
				//safe penalized edge values, and set them to original values
				distances2_copy.clear();
				for (int i=0; i<customers; i++)
				{
			    	distances2_copy.add((HashMap<Integer, Float>) nodes[i].distances_p.clone());
			    	nodes[i].distances_p = (HashMap<Integer, Float>) nodes[i].distances.clone();
				}

				//start local search moves from all routes that were changed during perturbation
				for (int i=0; i<routes_changed_during_perturbation.size(); i++)
				{
					int r = routes_changed_during_perturbation.get(i);
					for (int j=1; j<routes.get(r).length+1; j++)
					{
						starts_CROSS.add( new Edge(routes.get(r).nodes.get(j-1).index, routes.get(r).nodes.get(j).index) );
						if (j>1)
							starts_RC.add(routes.get(r).nodes.get(j-1).index);
					}
				}				
				for(int i=0; i< routes.size(); i++)
					routes_update_RC_variables.add(i);	    	
				penalties_off = true;
			}//reoptimisation	    
		
			//---------------------------------------------------------
			//EDGE PENALISATION
			//---------------------------------------------------------			 
			if (!penalties_off)
			{
				stopwatchOperators.start();
				//find the "worst" edge
				updateEdgeValues();
				Edge worstEdge = new Edge (-999999, -1, -1);
				for (int i=0; i<customers; i++)
				{
					if (edgeValues[i][0].costs > worstEdge.costs )
						worstEdge = edgeValues[i][0];
					if (edgeValues[i][1].costs > worstEdge.costs )
						worstEdge = edgeValues[i][1];
				}
		
				 if( worstEdge.node1 == -1 )
				 {
					 System.out.println("ERROR : No edge could be penalized"); 
					 break;
				 }
				 else
				 {
					 //penalize edge
					 int edgeNode1 = worstEdge.node1;
					 int edgeNode2 = worstEdge.node2;
					 int pen = getPenalties(edgeNode1, edgeNode2);
					 setPenalties(edgeNode1, edgeNode2, pen+1);
					 setDistance2(edgeNode1, edgeNode2, (float)(getDistance2(edgeNode1, edgeNode2) + penalty));
					 

					 //update the penalisation value of this edge
					 if (edgeNode1 < customers &&  edgeValues[edgeNode1][0].node2 == edgeNode2 ) 
						 edgeValues[edgeNode1][0].costs = edgeValues[edgeNode1][0].costs  * (pen+1) / (pen+2);
					 else if (edgeNode1 < customers)
						 edgeValues[edgeNode1][1].costs = edgeValues[edgeNode1][1].costs* (pen+1) / (pen+2);
					 if (edgeNode2 < customers && edgeValues[edgeNode2][0].node2 == edgeNode1) 
						 edgeValues[edgeNode2][0].costs = edgeValues[edgeNode2][0].costs* (pen+1) / (pen+2);
					 else if (edgeNode2 < customers)
						 edgeValues[edgeNode2][1].costs = edgeValues[edgeNode2][1].costs* (pen+1) / (pen+2);

					 //update detours and insertion costs
					 nodes[edgeNode1].detour_p += penalty;
					 nodes[edgeNode2].detour_p += penalty;
					 //update the values of the nearest neighbours
					 for (int nn=0; nn<nodes[edgeNode1].is_near_customer_of.size(); nn++)
						updateNearestNeighbour(nodes[edgeNode1].is_near_customer_of.get(nn), nodes[edgeNode1]);	
					 for (int nn=0; nn<nodes[edgeNode2].is_near_customer_of.size(); nn++)
						updateNearestNeighbour(nodes[edgeNode2].is_near_customer_of.get(nn), nodes[edgeNode2]);
		
					 //try to remove this edge in the next iteration with local search
					 if (edgeNode1<customers) starts_RC.add(edgeNode1);
					 if (edgeNode2<customers) starts_RC.add(edgeNode2);	
					 starts_CROSS.add( new Edge(edgeNode1, edgeNode2) );
					 penalties_given++;
					 penalizedEdges.add(new Edge(edgeNode1,edgeNode2));
			 }
			 stopwatchOperators.stop();
			 run_time_operators[3] += stopwatchOperators.elapsedTimeMillis();		
			} 
			
	    }//iterations
	    
		//---------------------------------------------------------
		//FINALIZE SEARCH
		//---------------------------------------------------------
	 // System.out.println("computation times : " + run_time_operators[0] + ", " + run_time_operators[1] + ", " +run_time_operators[2] + ", " +run_time_operators[3]);
	 // System.out.println(countPatternTotal + ", " + countPatternNotCurrent + ", " + countPatternImprove);
	  
	    resetSolution( best_routing, null, null);
	    //getNearestNeighbourStats();
	    //paint the final solution
	    ArrayList<ArrayList<Integer>> routesCopy = new ArrayList<ArrayList<Integer>>();
    	for (int i=0; i<routes.size(); i++)
    	{
    		routesCopy.add( new ArrayList<Integer>() );
			for (int j=0; j<=routes.get(i).length+1; j++)
				routesCopy.get(i).add( routes.get(i).nodes.get(j).index );
    	}

	}
	
	private void resetSolution(ArrayList<Route> targetState, float[][] penalizedDistance, int[][] penalties)
	{
		//reset the routes and customers
		 for (int r=0; r<routes.size(); r++)
		 {
		    routes.get(r).cost 	= targetState.get(r).cost;
		    routes.get(r).load 	= targetState.get(r).load;
		    routes.get(r).length= targetState.get(r).length;
		    routes.get(r).depot = targetState.get(r).depot;
		    routes.get(r).nodes.clear();
		    for (int n=0; n<targetState.get(r).length+2; n++)
		    {
		    	routes.get(r).nodes.add( targetState.get(r).nodes.get(n) );
		    	if (targetState.get(r).nodes.get(n).index<customers)
		    	{
			    	routes.get(r).nodes.get(n).routeIndex  = r;
			    	routes.get(r).nodes.get(n).routePosition = n;
			    	routes.get(r).nodes.get(n).depot = routes.get(r).depot;
		    	}
		    }
		 }				 
		 //remove previous penalization
		 /*if ( penalties != null)
		 {
			//reset penalization to a previous state
			for (int i=0; i<nodes.length; i++)
			{
				nodes[i].distances2 = penalizedDistance[i].clone();
				nodes[i].penalties 	= penalties[i].clone();
			}
		 }
		else*/
		{
			//remove all penalties			
			for (int i=0; i<nodes.length; i++)
			{			
				Iterator it = nodes[i].penalties.entrySet().iterator();
				while (it.hasNext()) 
				{
				       Map.Entry pair = (Map.Entry)it.next();
				       pair.setValue(0);
				}
				
				it = nodes[i].distances_p.entrySet().iterator();
				while (it.hasNext()) 
				{
				       Map.Entry pair = (Map.Entry)it.next();
				       pair.setValue(nodes[i].distances.get(pair.getKey()));
				}

			}
			
		}
		 for(int i=0; i< routes.size(); i++)
		 {
			routes_to_be_optimised.add(i);//to update neighbours
			routes_update_edge_values.add(i);//to update ejection chain memory
			routes_update_RC_variables.add(i);//to update edge values
		 }
		 update_neighbours();
		 routes_to_be_optimised.clear();
		 
		 //update all relevant variables after the reset
		// routeChanges			= 0;
		 //routeChanges2 			= 0;
		 //reoptimisation_count 	= 1;//for reOpt of entire solution
	}
	
	
	private static class MapUtil {
	    public static <K, V extends Comparable<? super V>> Map<K, V> 
	        sortByValue(Map<K, V> map) {
	        List<Map.Entry<K, V>> list = new LinkedList<Map.Entry<K, V>>(map.entrySet());
	        Collections.sort( list, new Comparator<Map.Entry<K, V>>() {
	            public int compare(Map.Entry<K, V> o1, Map.Entry<K, V> o2) {
	                return (o1.getValue()).compareTo( o2.getValue() );
	            }
	        });

	        Map<K, V> result = new LinkedHashMap<K, V>();
	        for (Map.Entry<K, V> entry : list) {
	            result.put(entry.getKey(), entry.getValue());
	        }
	        return result;
	    }
	}
	
	private void sort_insertions()
	{
		for (int i=0; i<customers; i++)
			nodes[i].insert_next_to = (HashMap<Integer, InsertionCost>) MapUtil.sortByValue(nodes[i].insert_next_to);
	}
	
	public boolean isFeasibleMTVRP(int maxVehicles, double maxCosts)
	{
		//1.) get a list of the cost values of all routes
		int activeRoutes = 0;
	    for (int r=0; r<routes.size(); r++)
			 if (routes.get(r).length > 0)
				 activeRoutes++;
	    ArrayList<Double> costList = new ArrayList<Double>();
	    int index = 0;
	    double totalCosts =0;
	    for (int r=0; r<routes.size(); r++)
			 if (routes.get(r).length > 0)
			 {
				 costList.add( routes.get(r).cost );
			 }
	    
	    ArrayList<Double>  costListOriginalFF = (ArrayList<Double>) costList.clone();
	    ArrayList<Double>  costListOriginal = (ArrayList<Double>) costList.clone();
	    Collections.sort(costList);
	    Collections.reverse(costList);
	    if (costList.get(0) > maxCosts)
	    	;//System.out.println("not possible");
	    
	    //System.out.println(routes.get(0).cost + " " + routes.get(1).cost);
	    
	    
	    //Heuristic for bin-packing: largest possible first FIRST FIT DECREASING
	    int countVehicles = 1;
	    double spaceCurrentVehicle = maxCosts;
	    while ( !costList.isEmpty() )
	    {
	    	int pos = 0;
	    	while ( costList.get(pos) > spaceCurrentVehicle )
	    	{
	    		pos++;
	    		if ( pos == costList.size())//need next vehicle
	    		{
	    			countVehicles++; spaceCurrentVehicle=maxCosts; pos = 0; break;
	    		}
	    	}
	    	spaceCurrentVehicle -= costList.get(pos);
	    	costList.remove(pos);
	    }
	    if (countVehicles <= maxVehicles)
	    	return true;
	    
	  //Heuristic for bin-packing: largest possible first FIRST FIT 
	    countVehicles = 1;
	    spaceCurrentVehicle = maxCosts;
	    while ( !costListOriginalFF.isEmpty() )
	    {
	    	int pos = 0;
	    	while ( costListOriginalFF.get(pos) > spaceCurrentVehicle )
	    	{
	    		pos++;
	    		if ( pos == costListOriginalFF.size())//need next vehicle
	    		{
	    			countVehicles++; spaceCurrentVehicle=maxCosts; pos = 0; break;
	    		}
	    	}
	    	spaceCurrentVehicle -= costListOriginalFF.get(pos);
	    	costListOriginalFF.remove(pos);
	    }
	    if (countVehicles <= maxVehicles)
	    	return true;
	    
	    //MAX REST
	    double[] spaceLeft = new double[maxVehicles];
	    Arrays.fill(spaceLeft, maxCosts);
	    while ( !costListOriginal.isEmpty() )
	    {
	    	//find bucket with most remaining capacity
	    	double maxSpace = 0;
	    	int smallestIndex = -1;
	    	for (int i=0; i<maxVehicles; i++)
	    		if (spaceLeft[i]>maxSpace)
	    		{
	    			maxSpace = spaceLeft[i]; smallestIndex = i;
	    		}
	    	if (costListOriginal.get(0) <= maxSpace)
	    	{
	    		spaceLeft[smallestIndex] -= costListOriginal.get(0);
	    		costListOriginal.remove(0);
	    		if ( costListOriginal.isEmpty() )
	    			return true;
	    	}	    		
	    	else
	    		break;
	    }
	 
	    	return false;
	}
	
	public void extractPattern(int minPatternSize, int maxPatternSize, boolean applyMining)
	{
		ArrayList<Long> patternList = new ArrayList<Long>();
		ArrayList<Integer> pattern = new ArrayList<Integer>();
		for (int patternSize = minPatternSize; patternSize <= maxPatternSize; patternSize++)
		for (Route r : routes)
		{
			if (r.length > 0)//ignore empty routes
				for (Node n : r.nodes)
				{
					if (r.length +2 >= n.routePosition + patternSize && r.length >= patternSize)//are there still sufficient nodes in the route to extract the pattern?
					//r.length indicates the number of customers in the route
					{
						//extract pattern
						pattern.clear();
						int currentSize = 0;
						while (currentSize < patternSize)
						{
							pattern.add(r.nodes.get(n.routePosition + currentSize).index);
							currentSize++;
						}
						//reverse the pattern, if that results in a smaller identifier
						if (pattern.get(0) < pattern.get(pattern.size() - 1))
							Collections.reverse(pattern);
						//compute the identifier
						long identifier = 0;
						for (int i = 0; i < patternSize; i++)
							identifier += (long)Math.pow((long)customers+1, (long)i) * (long)pattern.get(i);
						
						long tempID = identifier;
					
						//increase the count value of the identifier
						int count = minedPatterns.containsKey(identifier) ? minedPatterns.get(identifier) : 0;
						minedPatterns.put(identifier, count + 1);	
						patternList.add( identifier );
						
						if ((count == minMaxFrequency || mostFrequPattern.size() < numConsideredPattern) && (!mostFrequPattern.containsKey(identifier))){
							addHighFrequPattern(identifier, count+1);
						}
						else if (count > minMaxFrequency && mostFrequPattern.containsKey(identifier) )
							mostFrequPattern.put(identifier, count + 1);
					}
					else
						break;
						
				}
		}
		patternHistory.add( patternList );
		if (patternHistory.size() > 100 )// only store the pattern of the last X local minima
		{
			for (long id : patternHistory.get(0))
			{
				int count = minedPatterns.get(id);
				if (count == 1){
					minedPatterns.remove(id);
					mostFrequPattern.remove(id);
				}
				else{
					minedPatterns.put(id, count - 1);
					mostFrequPattern.put(id, count - 1);
				}
			}
			patternHistory.remove( 0 );
		}
		
		//ignore the last 10
		//fill once the history is complete
		/*if (patternHistory.size() == 19 )
		{
			for (int i =0; i<10;i++)
			for (long id : patternHistory.get(i))
			{
				int count = minedPatterns.containsKey(id) ? minedPatterns.get(id) : 0;
				minedPatterns.put(id, count + 1);	
			}
		}
		//update
		if (patternHistory.size() == 20 )
		{
			for (long id : patternHistory.get(9))
			{
				int count = minedPatterns.containsKey(id) ? minedPatterns.get(id) : 0;
				minedPatterns.put(id, count + 1);	
			}
		}*/
		applyMining( patternList );
		//System.out.println(minedPatterns.size());
	}
	
	
	public void applyMining ( ArrayList<Long> patternList )
	{
		for (Map.Entry<Long, Integer> entry : mostFrequPattern.entrySet()) 
		{
			countPatternTotal++;
		    long identifier = entry.getKey();
		    if (!patternList.contains(identifier))
		    {
		    	//countPatternNotCurrent++;
		    	//int frequency = (int) Math.round(10.0 *(double)(int) entry.getValue() / (double)maxCollection );
		    	//averageStats[ frequency ]++;		    	
		    	//averageMoveChanges[getChanges ( identifier )][frequency]++;
		    	//averageRouteChanges[ getRouteChanges( identifier)][frequency]++;
		    	evaluatePattern2(identifier);
		    }    
		}
		/*statsCollected++;
		if (statsCollected==Math.round( 300/maxCollection ))
		{
			int sum = 0;
			for (int i=10; i>-1; i--)
			{
				sum +=averageStats[i];
				System.out.println( Math.round((double)sum / (double)statsCollected));
			}
			sum = 0;
			for (int i=10; i>-1; i--)
			{
			for (int j=9; j>-1; j--)
			{
				sum +=averageMoveChanges[j][i];
				System.out.println( j + "changes " + Math.round((double)averageRouteChanges[j][i] / (double)statsCollected));
			}
			System.out.println("-------------------"); 
			}
			
			System.out.println("---------------------");
			sum = 0;
			for (int i=10; i>-1; i--)
			{
				sum +=averageNumNodes[i];
				System.out.println( Math.round((double)sum / (double)statsCollected));
			}
		}
		patternCollected = 0;
		minedPatterns.clear();
		System.out.println("aha");	*/
	}
	
	public int getChanges(long identifier)
	{
		ArrayList<Integer> sequence = new ArrayList<Integer>();
		for (int i=maxPatternSize-1; i>-1; i--)//obtain the nodes from the identifier
	    {
	    	int nextNode = (int) Math.floor( identifier / Math.pow(customers+1, i) );
	    	identifier -= nextNode * Math.pow(customers+1, i);
	    	sequence.add( nextNode );
	    }
		int countEdgeBreaks = 0;
		ArrayList<ArrayList<Integer>> newNeighbours = new ArrayList<ArrayList<Integer>>();
		for (int i=0; i<maxPatternSize; i++)
		{
			if (i>0)
				if ( nodes[sequence.get(i)].neighbours[0].index != sequence.get(i-1) && nodes[sequence.get(i)].neighbours[1].index != sequence.get(i-1) )
					countEdgeBreaks++;
			if (i<maxPatternSize-1)
				if ( nodes[sequence.get(i)].neighbours[0].index != sequence.get(i+1) && nodes[sequence.get(i)].neighbours[1].index != sequence.get(i+1) )
					countEdgeBreaks++;			
		}
		return countEdgeBreaks;
		
	}
	
	public int getRouteChanges(long identifier)
	{
		ArrayList<Integer> involvedRoutes = new ArrayList<Integer>();
		for (int i=maxPatternSize-1; i>-1; i--)//obtain the nodes from the identifier
	    {
	    	int nextNode = (int) Math.floor( identifier / Math.pow(customers+1, i) );
	    	identifier -= nextNode * Math.pow(customers+1, i);
	    	if (!involvedRoutes.contains( nodes[nextNode].routeIndex) )
	    		involvedRoutes.add( nodes[nextNode].routeIndex );
	    }
		return involvedRoutes.size();
		
	}
	
	
	public void evaluatePattern2(long identifier)
	{
		//------------------------------------------
		//extract the pattern from the identifier
		//------------------------------------------
		ArrayList<Integer> pattern = new ArrayList<Integer>();
		ArrayList<Route> involvedRoutes = new ArrayList<Route>();
		int patternSize = 1;
		while (identifier > Math.pow(customers+1, patternSize))
			patternSize++;
		
		for (int i=patternSize-1; i>-1; i--)//obtain the nodes from the identifier
	    {
	    	int nextNode = (int) Math.floor( identifier / Math.pow(customers+1, i) );
	    	if (nextNode == customers)//is depot node
				return;//pattern with the depot node are currently not supported
	    	identifier -= nextNode * Math.pow(customers+1, i);
	    	pattern.add( nextNode );
	    	int r = nodes[nextNode].routeIndex;
	    	if (!involvedRoutes.contains(routes.get(r)))
	    		involvedRoutes.add(routes.get(r));
	    }
		
		//------------------------------------------
		//split the involved routes into fragments that need to be reconnected
		//------------------------------------------		
		ArrayList<Integer> reconnectionPoints = new ArrayList<Integer>();//nodes that need a new neighbour
		ArrayList<EdgeSequence> routeFragments = new ArrayList<EdgeSequence>();//fragments
		double sum_broken_Edges = 0;//the costs that are saved due to removed edges
		int depot_counter = 0;//whenever the depot node is reached, this is increases to have an unique id for the depot
		for (Route r : involvedRoutes)
		{
			ArrayList<Integer> fragment = new ArrayList<Integer>();
			boolean continue_fragment = false;
			int sum_demand = 0;
			//the start and end of a route (the depot node) is never part of the pattern and always part of a fragment
			for (Node n : r.nodes)
			{
				if ( !pattern.contains(n.index))//the current node is not part of the pattern
				{
					if ( continue_fragment )//continue the current fragment
					{
						fragment.add(n.index);
						sum_demand += n.demand;
					}						
					else//start a new fragment
					{
						fragment.add(n.index < customers ? n.index : n.index + depot_counter++);
						sum_demand += n.demand;
						continue_fragment = true;						
					}
				}
				else //the current node is part of the pattern
				{
					//case 1: the previous node is part of a fragment, the current one of a pattern: fragment ends, pattern starts 
					if (continue_fragment)
					{
						//store the fragment up to this node
						routeFragments.add(new EdgeSequence(sum_demand, fragment.get(0), fragment.get(fragment.size()-1), (ArrayList<Integer>) fragment.clone()));						
						reconnectionPoints.add(fragment.get(fragment.size()-1));
						//clear the fragment
						fragment.clear();
						sum_demand = 0;
						continue_fragment = false;
						//break the edge to the previous node
						sum_broken_Edges += getDistance(n.index, n.neighbours[0].index);
					}
					//case 2: the previous node is also part of the pattern, but not an immediate neighbour in the pattern, i.e., there is no edge between them
					else if ( Math.abs( pattern.indexOf( n.index ) - pattern.indexOf( n.neighbours[0].index ) ) != 1 )
					{
						sum_broken_Edges += getDistance(n.index, n.neighbours[0].index);
					}
					//is the successor part of the pattern? If not, the edge is broken
					if ( !pattern.contains(n.neighbours[1].index) )//is the edge to the succesor broken?
					{
						sum_broken_Edges += getDistance(n.index, n.neighbours[1].index);
						reconnectionPoints.add(n.neighbours[1].index < customers ? n.neighbours[1].index : n.neighbours[1].index + depot_counter);
					}
				}						
			}
			//finish the last fragment
			int end_of_fragment = fragment.get(fragment.size()-1) < customers ? fragment.get(fragment.size()-1) : customers + depot_counter++;
			routeFragments.add(new EdgeSequence(sum_demand, fragment.get(0), end_of_fragment, (ArrayList<Integer>) fragment.clone()));
		}
		
		//recode and store pattern as fragment 
		int demandPattern = 0;
		for (int i : pattern)
			demandPattern += nodes[i].demand;
		routeFragments.add(new EdgeSequence(demandPattern, pattern.get(0), pattern.get(pattern.size()-1), pattern));
		reconnectionPoints.add(pattern.get(0)); reconnectionPoints.add(pattern.get(pattern.size()-1));
				
		//map reconnectionPoints to sequences
		HashMap<Integer, EdgeSequence>  mapping = new HashMap<Integer, EdgeSequence>();
		for (int r:reconnectionPoints)
			for (EdgeSequence s: routeFragments)
				if (s.startNode == r || s.endNode == r)
				{
					mapping.put(r, s);
					break;
			   }
		
		//detect new edges that occur in the pattern and deduct their costs from the savings
		for (int i=0; i<pattern.size()-1; i++)
			if ( nodes[pattern.get(i)].neighbours[0].index != pattern.get(i+1) && nodes[pattern.get(i)].neighbours[1].index != pattern.get(i+1) )
				{sum_broken_Edges -= this.getDistance(pattern.get(i), pattern.get(i+1));}
	
				
		//get all combination of pairs of connection points
		List<Integer> set = new LinkedList<Integer>(reconnectionPoints);
		ArrayList<List<List<Integer>>> all_combinations = new ArrayList<List<List<Integer>>>();		
		
		if (set.size()<10)//TODO limit number of reconnection parts
			getAllPairs(set, new ArrayList<List<Integer>>(), all_combinations, sum_broken_Edges);
		
		if (true)
		{
		//go through each possible combination and check for feasibility
		//1) no cycles in routes, each route starts and ends at a depot
		//2) capacity limit is not violated
		//3) the sum of the costs of the new edges minus the sum of the costs of the removed edges is negative
		for (List<List<Integer>> candidate : all_combinations)
		{
			double sum_new_edges = 0;//the costs of the new edges
			boolean hasCycles = false;
			//check whether there is a positive gain in implementing the pattern and no sub-tours
			//TODO: this can be done in a sequential search fashion!
			for (List<Integer> pair : candidate)
			{
			    if (mapping.get(pair.get(0)) == mapping.get(pair.get(1)))//no connection between reconnection points of the same fragment
			    	hasCycles = true;
			    if (mapping.get(pair.get(0)).startNode < customers &&  mapping.get(pair.get(0)).endNode < customers && mapping.get(pair.get(1)).startNode < customers &&  mapping.get(pair.get(1)).endNode < customers)//no connection between two subsequences TODO: allow them
			    	hasCycles = true;//check this condition!
			    if (mapping.get(pair.get(0)).accCapacity + mapping.get(pair.get(1)).accCapacity > this.capacityLimitVehicle)
			    	hasCycles = true;	
			    sum_new_edges += getDistance(Math.min(customers, pair.get(0)), Math.min(customers, pair.get(1)));
			 }
			
			//construct the routes and do the capacity check
			//TODO the code below depends on how tours are represented, the code can be improved
			 if ( sum_broken_Edges > sum_new_edges && !hasCycles)
			 {			    
				 //the new routes
				 ArrayList<ArrayList<Integer>> constructedRoutes = new ArrayList<ArrayList<Integer>>();
				 ArrayList<Integer> capacity_routes = new ArrayList<Integer>();
				 //reconnect all fragments
				 for (List<Integer> pair : candidate)
				 {
					boolean connection_extended = false;
					//check whether one of the fragments has already been added to a route
					for (ArrayList<Integer> r : constructedRoutes)
					{
						//case 1: attach second fragment to the front
						if ( r.get(0) == (int)pair.get(0) )
						{
							if ( mapping.get(pair.get(1)).startNode == (int)pair.get(1) )
								Collections.reverse( mapping.get(pair.get(1)).nodes );
							r.addAll(0, mapping.get(pair.get(1)).nodes);
							connection_extended = true;
						}
						//case 2: attach first fragment to the front
						else if ( r.get(0) == (int)pair.get(1) )
						{
							if ( mapping.get(pair.get(0)).startNode == (int)pair.get(0) )
								Collections.reverse( mapping.get(pair.get(0)).nodes );
							r.addAll(0, mapping.get(pair.get(0)).nodes);
							connection_extended = true;
						}
						//case 3: attach first fragment to the back
						else if ( r.get(r.size()-1) == (int)pair.get(0) )
						{
							if ( mapping.get(pair.get(1)).endNode == (int)pair.get(1) )
								Collections.reverse( mapping.get(pair.get(1)).nodes );
							r.addAll(mapping.get(pair.get(1)).nodes);
							connection_extended = true;
						}
						//case 4: attach second fragment to the back
						else if ( r.get(r.size()-1) == (int)pair.get(1) )
						{
							if ( mapping.get(pair.get(0)).endNode == (int)pair.get(0) )
							{
								Collections.reverse( mapping.get(pair.get(0)).nodes );
							}
							r.addAll(mapping.get(pair.get(0)).nodes);
							connection_extended = true;
						}
					}
					//case 5: connect both fragments and store the result as a partial route
					if (!connection_extended)
					{
						if (mapping.get(pair.get(0)).startNode == (int)pair.get(0))
							Collections.reverse( mapping.get(pair.get(0)).nodes );
						if (mapping.get(pair.get(1)).endNode == (int)pair.get(1))
							Collections.reverse( mapping.get(pair.get(1)).nodes );
						ArrayList<Integer> newRoute = new ArrayList<Integer>();
						newRoute.addAll(mapping.get(pair.get(0)).nodes);
						newRoute.addAll(mapping.get(pair.get(1)).nodes);
						constructedRoutes.add( newRoute );
					}
					
					/*System.out.println(pair.get(0) + " - "  + pair.get(1) + " " );
					for (ArrayList<Integer> r : new_routes_int)
					 {
						System.out.println( r.get(r.size()-1) );
						 for (int n: r)
							 System.out.print(n+"-");
						 System.out.println();
					 }
					 System.out.println("----");*/
				 }
				 
				 /*for (EdgeSequence n : routeFragments)
					{
						System.out.println ( "(" +  n.startNode+","+ n.endNode+ ") " );
					}*/
				 
				 //count the number of new edges
				 
				/* for (Integer n : pattern)
						System.out.print (nodes[n].routeIndex+"("+nodes[n].routePosition+") - ");
					System.out.println();
					for (EdgeSequence n : sequences)
					{
						Node temp = nodes[customers]; 
						if (n.startNode < customers)
							temp = nodes[n.startNode]; 
						System.out.println ( "(" +  n.startNode+","+ n.endNode+ ") "  + "  " + n.accCapacity);
					}
					for (Integer i : reconnectionPoints)
						System.out.println (i);
					for (List<Integer> pair : candidate)
						System.out.println(pair);
					System.out.println("KKKKKKKKKKKKKKKKKKK  " + sum_broken_Edges);*/
				 
				 //construct routes
				 ArrayList<Route> newRoutes = new ArrayList<Route>();			 
				 for (ArrayList<Integer> r : constructedRoutes)
				 {
					 Route route = new Route();
					 for (int i=0; i<r.size(); i++)
					 {
						 Node n = nodes[Math.min(r.get(i), customers)];
						 route.nodes.add(n);
						 route.load += n.demand;
					 }
					 route.length = route.nodes.size() - 2;
					 if ( route.load > capacityLimitVehicle )
						 return;
					 newRoutes.add( route );
				 }
				 
				 int countNewRoutes = 0;
				 for (Route r : newRoutes)
				 {
					 for (int i=0; i<r.nodes.size(); i++)
					 {
						 Node n = r.nodes.get(i);
						 if (n.index < customers)
						 {
							 n.routePosition = i;
						 }
					 }
					 countNewRoutes++;
				 }
				 
				 /*for (ArrayList<Integer> r : constructedRoutes)
					{
						for (Integer n: r)
							System.out.print(n + " - ");
						System.out.println();
					}
					System.out.println("---");
					for (Route r : involvedRoutes)
					{
						for (Node n: r.nodes)
							System.out.print(n.index + " - ");
						System.out.println();
					}
					System.out.println("---");
					*/
				
					
				 int newEdges = 0;
				 for (Route r : newRoutes)
				 {
					 for (int i=1; i<r.nodes.size()-1; i++)
					 {
						 if (r.nodes.get(i).neighbours[0] != r.nodes.get(i-1) && r.nodes.get(i).neighbours[1] != r.nodes.get(i-1))
							 newEdges++;
					 }
					 //last customer to depot
					 if (r.nodes.get(r.length).neighbours[0] != r.nodes.get(r.length+1) && r.nodes.get(r.length).neighbours[1] != r.nodes.get(r.length+1))
						 newEdges++;
				 }
				 
				 if ( involvedRoutes.size()==1)
					;//System.out.println("intra " + newEdges);
				 else
					;// System.out.println(involvedRoutes.si);
				// System.out.println(newEdges);
					//System.out.println("-------------------");
					
				 averageStatsMoves[newEdges]++;
				 

				 
				 double values_before = this.computeSolutionCosts();
				 //remove old routes
				 for (Route r:involvedRoutes)
					 routes.remove(r);
				 //add new ones
				 routes.addAll( newRoutes );
				 for (int i = 0; i< routes.size(); i++)
					 for (int j=1; j<routes.get(i).nodes.size()-1; j++)
					 {
						 routes.get(i).nodes.get(j).routeIndex = i;
						 routes_to_be_optimised.add(i);
					 }
				 update_neighbours();
				 
				 double values_after = this.computeSolutionCosts();
				 this.checkSolution();
				 if (values_before - values_after != sum_broken_Edges - sum_new_edges)
					 System.out.println("wrong costs computed");
				 //System.out.println( values_before - values_after + ", " + sum_broken_Edges + " - " + sum_new_edges);
				 //System.out.println( "Moved found: " + involvedRoutes.size());
				 averageStats[pattern.size()]++;
				 					
					return;
			 }
		}	
		}
	}
	
	
	//private static void getAllPairs(Set<Integer> set,
	private void getAllPairs(List<Integer> set,
	        List<List<Integer>> currentResults,
	        List<List<List<Integer>>> results,
			double saving)	
	    {
	        if (set.size() < 2)
	        {
	            results.add(new ArrayList<List<Integer>>(currentResults));
	            return;
	        }
	        List<Integer> list = new ArrayList<Integer>(set);
	        Integer first = list.remove(0);
	        for (int i=0; i<list.size(); i++)
	        {
	            Integer second = list.get(i);
	            double connectionCost = getDistance(Math.min(customers, first), Math.min(customers, second));
	            //Set<Integer> nextSet = new LinkedHashSet<Integer>(list);
	            //nextSet.remove(second);
	            if (saving > connectionCost)
	            {
		            List<Integer> nextSet = new LinkedList<Integer>(list);
		            nextSet.remove(second);
	
		            List<Integer> pair = Arrays.asList(first, second);
		            currentResults.add(pair);
		            getAllPairs(nextSet, currentResults, results, saving - connectionCost);
		            currentResults.remove(pair);
	            }
	        }
	    }
	
	
	
	
	
	public void getNearestNeighbourStats(){
		ArrayList<Integer> rankList = new ArrayList<Integer>();

	for(int i=0; i<customers; i++)
	{
		ArrayList<Edge> sortedDists = new ArrayList<Edge>();
		for(int j=0; j<customers; j++)
		{
			if (j != i)//distance to other customer
				sortedDists.add( new Edge(this.getDistance(i, j),j,-1) );
		}	
		Collections.sort(sortedDists);
		Collections.reverse(sortedDists);
		int countNeighours = 0;
		if (nodes[i].neighbours[0].index>=customers)
			countNeighours++;
		if (nodes[i].neighbours[1].index>=customers )
			countNeighours++;
		int rank = 0;
		while ( countNeighours < 2)
		{
			if (nodes[i].neighbours[0].index==sortedDists.get(rank).node1 || nodes[i].neighbours[1].index==sortedDists.get(rank).node1 )
				countNeighours++;
			rank++;	
		}
		//int degree1 = -1; if (nodes[i].neighbours[0].index < customers) degree1 = nodes[i].distanceRank[nodes[i].neighbours[0].index];
		//int degree2 = -1; if (nodes[i].neighbours[1].index < customers) degree2 = nodes[i].distanceRank[nodes[i].neighbours[1].index];
		//int maxD	= Math.max(degree1, degree2);
		rankList.add(rank);
	}
	int ninetyFive = (int)( Math.ceil(customers * 0.95) ) - 1;
	int ninetyNine = (int)( Math.ceil(customers * 0.99) ) - 1;
	Collections.sort(rankList);
	System.out.println(customers + ": " + rankList.get(ninetyFive) + "   " + rankList.get(ninetyNine));
	
	
	//LARGE SCALE TESTING
	/*		int max = 0;
			ArrayList<Integer> top1= new ArrayList<Integer>();
			for (int i=0; i<customers; i++)
			{
				int degree1 = -1; if (nodes[i].neighbours[0].index < customers) degree1 = nodes[i].distanceRank[nodes[i].neighbours[0].index];
				int degree2 = -1; if (nodes[i].neighbours[1].index < customers) degree2 = nodes[i].distanceRank[nodes[i].neighbours[1].index];
				int maxD	= Math.max(degree1, degree2);
				//top1.size()<10)
					top1.add(maxD);
				//else if (top1.get(0) < maxD)
					//{top1.remove(0); top1.add(maxD);}
				//if (nodes[i].neighbours[0].index < customers &&  degree1> max)
					//max = nodes[i].distanceRank[nodes[i].neighbours[0].index];
				//if (nodes[i].neighbours[1].index < customers && nodes[i].distanceRank[nodes[i].neighbours[1].index] > max)
					//max = nodes[i].distanceRank[nodes[i].neighbours[1].index];
				
			}
			Collections.sort(top1);
			System.out.println(top1.get(0)+ "; " + top1.get(7) + "; " + top1.get(9));*/
	}
	
	public void addHighFrequPattern(long identifier, int count){
		
		if (mostFrequPattern.size() < numConsideredPattern){			
				//insert
				mostFrequPattern.put(identifier, count);
		}
		else{
			long smallestIdentifier = 0; 
			int smallestValue = Integer.MAX_VALUE;
			int secondSmallest = Integer.MAX_VALUE;
			Iterator it = mostFrequPattern.entrySet().iterator();
		    while (it.hasNext()) {
		        Map.Entry pair = (Map.Entry)it.next();
		        if ((int)pair.getValue() < smallestValue){
		        	smallestValue = (int)pair.getValue();
		        	smallestIdentifier = (long)pair.getKey();
		        }
		        else if ((int)pair.getValue() < secondSmallest){
		        	secondSmallest = (int)pair.getValue();
		        }
		    }
		    mostFrequPattern.remove(smallestIdentifier);
		    mostFrequPattern.put(identifier, count);
		    minMaxFrequency = secondSmallest; if (minMaxFrequency == Integer.MAX_VALUE ) minMaxFrequency = 0;    
		}
	}
	
	
	public MDVRPModel(JSONObject jsonObj) {
	
	    JSONArray customerInfo = (JSONArray) jsonObj.get("customers");     
	    Iterator<Object> iterator = customerInfo.iterator();
	    ArrayList<Node> nodes = new ArrayList<Node>();
	    int nodeIndex = 0;
	    while (iterator.hasNext()) {
	    	JSONObject nextEntry = (JSONObject) iterator.next();
	    	Node node = new Node(nodeIndex);
	    	node.demand = ((Long) nextEntry.get("demand")).intValue();
	    	node.demandPerProduct = new int[]{node.demand};
	    	node.x = ((Long) nextEntry.get("x")).intValue(); 
	    	node.y = ((Long) nextEntry.get("y")).intValue();
	    	nodes.add(node);
	    	nodeIndex++;
	    }
	
	    this.customers = nodes.size();
	    this.vertices = customers;
	    this.depots = 0;
	    this.capacityLimitVehicle = ((Long) jsonObj.get("vehicle_capacity")).intValue();
		this.cost_per_Route = ((Long) jsonObj.get("vehicle_costs")).doubleValue();
		this.nodes 	= new Node[customers];
		for (int i=0; i<customers; i++)
			this.nodes[i] = nodes.get(i);
		
		perturbation_length = 30;
		inventoryOff = false;
		distances_rounded = false;
		pruneLimit = Math.min(30, customers-1);
		
		computeAndSetDistanceMatrix();	
		
		this.products  		= 1;
		//lists
		routes 						= new ArrayList<Route>();
		routes_to_be_optimised 		= new ArrayList<Integer>();
		routes_changed_during_perturbation = new ArrayList<Integer>();
		starts_RC 					= new ArrayList<Integer>();
		routes_update_RC_variables 	= new ArrayList<Integer>();
		routes_update_edge_values 	= new ArrayList<Integer>();
		starts_CROSS     			= new ArrayList<Edge>();
		distance_to_nearest_depot  	=  new Edge[customers];
		this.maxTourLength = Double.MAX_VALUE;
	}
	
	
	private void computeAndSetDistanceMatrix() {

		int maxCWentries = customers<=300 ? customers: 100;
		for(int i=0; i<customers; i++) {
			ArrayList<Edge> sortedDists = new ArrayList<Edge>();
			for(int j=0; j<customers ; j++) {
				float dist = (float) Math.ceil(100.0 * Math.sqrt( Math.pow(nodes[i].x - nodes[j].x , 2) + Math.pow(nodes[i].y - nodes[j].y , 2)));		
				//float dist = (float) Math.sqrt( Math.pow(nodes[i].x - nodes[j].x , 2) + Math.pow(nodes[i].y - nodes[j].y , 2));	
				if (i != j)
					sortedDists.add( new Edge(dist,j,-1) );
			}	
			Collections.sort(sortedDists);
			Collections.reverse(sortedDists);

			for (int j=0; j<Math.min(sortedDists.size(), maxCWentries); j++)
				nodes[i].distances_for_CW.put(sortedDists.get(j).node1, (float) sortedDists.get(j).costs);
					
			for (int j=0; j<Math.min(sortedDists.size(), 1000); j++) {				
				nodes[i].distances.put(sortedDists.get(j).node1, (float) sortedDists.get(j).costs);
				nodes[i].distances_p.put(sortedDists.get(j).node1, (float) sortedDists.get(j).costs);
				nodes[i].penalties.put(sortedDists.get(j).node1, 0);
			}
			for (int j=0; j<pruneLimit; j++) {
				nodes[i].insert_next_to.put(sortedDists.get(j).node1, new InsertionCost(0,0,-1,-1));
				nodes[sortedDists.get(j).node1].is_near_customer_of.add(nodes[i]);//TODO where is this used?
				nodes[i].near_customers.add(sortedDists.get(j).node1);//TODO where is this used?
			}
			for (int j=0; j<Math.min(sortedDists.size(), 50); j++)
				nodes[i].customers_considered_as_neigbours.add(sortedDists.get(j).node1);
		}
    }
	
	
	public void setDepots(ArrayList<Node> depotConfiguration) {
		//reset model
		routes.clear();//TODO more resets?
		
		depots = depotConfiguration.size();
		vertices = customers + depots;
		
		//resize array containing all nodes
		nodes = Arrays.copyOf(nodes, vertices);
		depotInventory = new int[depots][1];
		for (int depotIndex=0; depotIndex<depots; depotIndex++) {
			Node node = depotConfiguration.get(depotIndex).copy();
			node.index = customers + depotIndex;
			depotInventory[depotIndex][0] = (int)node.demand;
			nodes[node.index] = node;
		}		
		
		//compute distances
		//compute the distances (only from customers)
		for (int i=0; i<customers; i++)
		{
			for(int j=customers; j<vertices; j++)
			{
				float dist = (float) Math.ceil(100.0 * Math.sqrt( Math.pow( nodes[i].x - nodes[j].x , 2) + Math.pow( nodes[i].y - nodes[j].y , 2) ) );
				//float dist = (float)  Math.sqrt( Math.pow( nodes[i].x - nodes[j].x , 2) + Math.pow( nodes[i].y - nodes[j].y , 2));
				
				//System.out.print(dist + " ");
				nodes[i].distances.put(j, (float) dist);
				nodes[i].distances_p.put(j, (float) dist);
				nodes[i].penalties.put(j, 0);
			}
			//System.out.println();
		}
		
	}
	
	
	public double constructStartingSolution() {
		
		/*depotInventory = new int[depots][1];						
		for (int i=0; i<depots; i++) {
			depotInventory[i][0] = (int)nodes[customers+i].demand;
		}*/
		ClarkWright cw= new ClarkWright();
		cw.clarkWright_MPMDVRPI( this, null, 999, false );
		return this.computeSolutionCosts();
	}

	
}//CLASS
