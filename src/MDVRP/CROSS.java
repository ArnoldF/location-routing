package MDVRP;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;


//demandPerProduct out-commented!
public class CROSS {

	 MDVRPModel model;
	
	private  class CROSS_start
	{		
		double savings;
		int directionRoute2;
		double dist_n2_break;
		double dist_n2_break_T;
		double dist_n1B_n2B;
		int setEmpty;
		
		public CROSS_start (double s, int d2, double dist1, double dist2, double dist3, int sE)
		{
			savings 		= s;
			directionRoute2	= d2;
			dist_n2_break	= dist1;
			dist_n2_break_T	= dist2;
			dist_n1B_n2B	= dist3;
			setEmpty		= sE;
		}
	}
	
	private  class CROSS_option implements Comparable<CROSS_option>
	{
		double savings;
		int route1;
		int route2;
		int set1_start;
		int set1_end;
		int set2_start;
		int set2_end;
		int capacityChangesRoute1;
		int capacityChangesRoute2;
		double costChangesRoute1;
		double costChangesRoute2;
		int directionSet1;
		int directionSet2;
		int emptySet_pos;
		
		public CROSS_option (double s, int r1, int r2, int set1s, int set1e, int set2s, int set2e, int emptySet)
		{
			savings 	= s;
			route1		= r1;
			route2		= r2;
			set1_start	= set1s;
			set1_end	= set1e;
			set2_start	= set2s;
			set2_end	= set2e;
			emptySet 	= emptySet_pos;
		}
		
		@Override
		public int compareTo(CROSS_option opt) 
		{
			if(opt.savings < this.savings)
				return -1;
			else if (opt.savings > this.savings)
				return 1;
			else
				return 0;
		}
	}
	
	
	
	@SuppressWarnings("unchecked")
	public  void cross(MDVRPModel MyModel, ArrayList<Edge> crossEdges)
	{
		model = MyModel;		
		ArrayList<CROSS_option> cross_options		=  new ArrayList<CROSS_option>();
		
		for (int edge = 0; edge < crossEdges.size(); edge++)
			for (int edgeNode = 0; edgeNode<2; edgeNode++)
			{
				//find starts that remove the edge, and connect node1 to a node in another route
				Node node1 = null;
				if (edgeNode == 0)
					node1 = model.nodes[crossEdges.get(edge).node1];	
				else
					node1 = model.nodes[crossEdges.get(edge).node2];
										
				if (node1.index < model.customers)
				{
				int node1_pos = node1.routePosition;
				int r1 = node1.routeIndex;
				Route route1 = model.routes.get(r1);
				
				//find near nodes of node1 (nn) that are in a different route
				for (int nn = 0; nn < model.pruneLimit; nn++)//node1.nearest.size()
					if ( model.nodes[node1.near_customers.get(nn)].routeIndex != node1.routeIndex)
					{			
						//------------------------------------------------------------------------
						//FIND CROSS-exchange starts with a positive saving
						//------------------------------------------------------------------------
						//1) safe variables for faster access
						//node2 is a node close to node 1 in a different route route2
						ArrayList<CROSS_start> possibleStarts = new ArrayList<CROSS_start>();
						Node node2 = model.nodes[node1.near_customers.get(nn)];
						Route route2 = model.routes.get(node2.routeIndex);//TODO what if node2 is depotNode?
						int node2_pos = node2.routePosition;
						int r2 = node2.routeIndex;
						int depot1 = route1.depot;
						int depot2 = route2.depot;
						
						//get the adjacent node of node1 of the edge that is to be removed
						Node node1_neighbour = null;
						if (edgeNode == 0)
							node1_neighbour = node1.neighbours[1];
						else
							node1_neighbour = node1.neighbours[0];
							
						//get the adjacent nodes of node2, if they are not depot nodes
						Node node2_succ = null; if (node2_pos<route2.length+1) node2_succ = node2.neighbours[1];
						Node node2_pre 	= null; if (node2_pos>0) 			   node2_pre = node2.neighbours[0];
						
						//2) get the distances of involved edges of the start
						//edges that are removed
						double dist_n1_break_p 		= model.getDistance2(node1.index, node1_neighbour.index);
						double dist_n1_break 		= model.getDistance(node1.index, node1_neighbour.index);
						double dist_n2_break_pre_p 	= -999999;
						double dist_n2_break_pre 	= model.getDistance(node2.index, node2_pre.index);
						if(node2_pos>0)
							dist_n2_break_pre_p 	= model.getDistance2(node2.index, node2_pre.index);
						double dist_n2_break_suc_p 	= -999999;
						double dist_n2_break_suc 	= model.getDistance(node2.index, node2_succ.index);
						if(node2_pos<route2.length+1)
							dist_n2_break_suc_p 		= model.getDistance2(node2.index, node2_succ.index);
						
						//edges that are added
						double dist_n1_n2 		= model.getDistance(node1.index, node2.index);
						double dist_n1N_n2N_pre	= 999999;
						//it is not allowed to build an edge between two depots
						if(node2_pre!=null && !(node2_pre.index>=model.customers && node1_neighbour.index>=model.customers))
							dist_n1N_n2N_pre = model.getDistance(node1_neighbour.index, node2_pre.index);
						double dist_n1N_n2N_suc	= 999999;
						//it is not allowed to build an edge between two depots
						if(node2_succ!=null && !(node2_succ.index>=model.customers && node1_neighbour.index>=model.customers))
							dist_n1N_n2N_suc = model.getDistance(node1_neighbour.index, node2_succ.index);					

						//3) find starts with a positive gain
						//option 1) remove edge between node2 and its predecessor
						double savingStart = dist_n1_break_p + dist_n2_break_pre_p - dist_n1_n2 - dist_n1N_n2N_pre;
						if ( savingStart >= 0 )
							possibleStarts.add(new CROSS_start(savingStart, -1, dist_n2_break_pre_p, dist_n2_break_pre, dist_n1N_n2N_pre, 0) );
						//option 2) remove edge between node2 and its successor
						savingStart = dist_n1_break_p + dist_n2_break_suc_p - dist_n1_n2 - dist_n1N_n2N_suc;
						if ( savingStart >= 0 )
							possibleStarts.add(new CROSS_start(savingStart, +1, dist_n2_break_suc_p, dist_n2_break_suc, dist_n1N_n2N_suc, 0) );							
											
						//option 3) a Crossover, where at least two nodes from route2 are integrated in rout1	
						//option 3a) the edge towards the predecessor is removed
						savingStart = dist_n1_break_p + dist_n2_break_pre_p - dist_n1_n2;
						if ( node2_pos < route2.length && savingStart - dist_n2_break_pre_p >= 0 && route1.load + node2.demand + node2_succ.demand <= model.capacityLimitVehicle )
							possibleStarts.add(new CROSS_start(savingStart, +1, dist_n2_break_pre_p, dist_n2_break_pre, -999999, 1) );
						//option 3b) the edge towards the successor is removed
						savingStart = dist_n1_break_p + dist_n2_break_suc_p - dist_n1_n2;
						if ( node2_pos > 1 && savingStart - dist_n2_break_pre_p >= 0 && route1.load + node2.demand + node2_pre.demand <= model.capacityLimitVehicle )
							possibleStarts.add(new CROSS_start(savingStart, -1, dist_n2_break_suc_p, dist_n2_break_suc, -999999, 1) );
						
						//option 4) a Crossover, where at least two nodes from route1 are integrated in route2
						savingStart = dist_n1_break_p + dist_n2_break_pre_p - dist_n1_n2;
						int dir = 1;
						if (edgeNode == 0)
							dir = -1;
						if ( node1_pos + dir <= route1.length && node1_pos + dir >= 1 && savingStart - dist_n2_break_pre_p > 0 && route2.load + node1.demand < model.capacityLimitVehicle )
							possibleStarts.add(new CROSS_start(savingStart, +1, dist_n2_break_pre_p, dist_n2_break_pre, -999999, 2) );
	
						//------------------------------------------------------------------------
						//FOR ALL STARTS, TRY TO FIND CROSS-exchange operations with a positive gain
						//------------------------------------------------------------------------
						while ( !possibleStarts.isEmpty() )
						{
							//CROSS-exchange moves							1
							if ( possibleStarts.get(0).setEmpty == 0 )
							{
								//two cases can be considered
								//cross_case 1:  node1 (and maybe more nodes) becomes part of route2
								//cross_case 2 : the neighbour of node1 becomes part of route2
								//the nodes in set1 go from route1 into route2
								//the nodes in set2 go from route2 into route1
								for (int cross_case = 1; cross_case <= 2; cross_case++ )
								{
									double saving_start =  possibleStarts.get(0).savings;
									int set1_end = -1;
									if ( cross_case == 1)
										set1_end = node1_pos;
									else
									{
										//the set cannot include a depot
										if (node1_neighbour.index >= model.customers)
											break;
										set1_end = node1_neighbour.routePosition;
									}

									//starting from a node in both routes, we extend the sets
									//the direction of extension has to be defined
									//-1 means 'backward' (we add the predecessor)
									//1  means 'forward' (we add the successor)
									int set2_direction = possibleStarts.get(0).directionRoute2;
									int set1_direction = 1;
									if (edgeNode == 0)
										set1_direction = -1;
									//the directions depend on which nodes change the routes, and the directions might have to be reversed
									if ( cross_case == 2)
									{
										set1_direction *= -1;
										set2_direction *= -1;
									}									
													
									//define the nodes with which the two sets start
									int set1_start = 0;
									int set2_start = 0;
									if ( cross_case == 1 )
									{
										set1_start = node1_pos;
										set2_start = node2_pos + set2_direction;
									}
									else
									{
										set1_start = node1_pos + set1_direction;
										set2_start = node2_pos;
									}
																		
									ArrayList<Integer> set1_nodes = new ArrayList<Integer>();
									set1_nodes.add( route1.nodes.get(set1_end).index );									
									int set1_cap = route1.nodes.get(set1_end).demand;
									int[] invNeeded1 = new int[model.products];
									for (int pr=0; pr < model.products; pr++)
										invNeeded1[pr] = route1.nodes.get(set1_end).demandPerProduct[pr];
									double set1_dist = 0;
										
									//extend set1 until the depot node is reached
									int lengthSet1 = 0;
									while (set1_end != 0 && set1_end != route1.length+1 && lengthSet1 < 500 )
									{
										lengthSet1++;
										int set2end = node2_pos + set2_direction;
										if ( cross_case == 2 )
											set2end = node2_pos;
										
										ArrayList<Integer> set2_nodes = new ArrayList<Integer>();
										set2_nodes.add( route2.nodes.get(set2end).index );
										int set2_cap =  route2.nodes.get(set2end).demand;
										int[] invNeeded2 = new int[model.products];
										for (int pr=0; pr<model.products; pr++)
											invNeeded2[pr] = route2.nodes.get(set2end).demandPerProduct[pr];
										double set2_dist = 0;									
										
										//extend set2 until the depot node is reached, keeping set1 fixed
										int lengthSet2 = 0;
										while (set2end != 0 && set2end != route2.length+1 && lengthSet2 < 500 &&
										route1.load - set1_cap + set2_cap <= model.capacityLimitVehicle)
										{
											lengthSet2++;
											//are the routes feasible if we exchange both sets?
											if ( route2.load - set2_cap + set1_cap <= model.capacityLimitVehicle
											&& (model.inventoryOff || depot1==depot2 || model.haveSufficientInventory(depot1, depot2, invNeeded1, invNeeded2)) 
											&& ( set1_nodes.size() < route1.length || set2_nodes.size() < route2.length ) )//do not exchange 2 complete routes
											{
												//get distance values of all edges that complete the move
												double dist_set1_end 			= model.getDistance(route1.nodes.get(set1_end).index, route2.nodes.get(set2end+set2_direction).index);
												double dist_set2_end 			= model.getDistance(route2.nodes.get(set2end).index, route1.nodes.get(set1_end+set1_direction).index);
												double dist_set1_end_break_p 	= model.getDistance2(route1.nodes.get(set1_end).index, route1.nodes.get(set1_end+set1_direction).index);
												double dist_set2_end_break_p 	= model.getDistance2(route2.nodes.get(set2end).index, route2.nodes.get(set2end+set2_direction).index);
												double dist_set1_end_break 		= model.getDistance(route1.nodes.get(set1_end).index, route1.nodes.get(set1_end+set1_direction).index);
												double dist_set2_end_break 		= model.getDistance(route2.nodes.get(set2end).index, route2.nodes.get(set2end+set2_direction).index);
												model.ev2++;
												//compute the changes to both routes in terms of length
												double changedDist_route1 = dist_set2_end + possibleStarts.get(0).dist_n1B_n2B - dist_set1_end_break - dist_n1_break - set1_dist + set2_dist;
												if ( cross_case == 2)
													changedDist_route1 = dist_set2_end + dist_n1_n2 - dist_set1_end_break - dist_n1_break - set1_dist + set2_dist;
												
												double changedDist_route2 = dist_set1_end + dist_n1_n2 - dist_set2_end_break - possibleStarts.get(0).dist_n2_break_T + set1_dist - set2_dist;
												if ( cross_case == 2)
													changedDist_route2 = dist_set1_end + possibleStarts.get(0).dist_n1B_n2B - dist_set2_end_break - possibleStarts.get(0).dist_n2_break_T + set1_dist - set2_dist;
												
												//compute the saving of the second cross
												double saving_end 		= dist_set1_end_break_p + dist_set2_end_break_p - dist_set2_end - dist_set1_end;
												double saving_CROSS		= saving_start + saving_end;
												
												//does the whole CROSS result in a saving, and is feasible in terms of route length?
												if ( saving_CROSS >= 0
														&& route1.cost + changedDist_route1 <= model.maxTourLength 
														&& route2.cost + changedDist_route2 <= model.maxTourLength 
														&& ( r1 > 2 || route1.cost + changedDist_route1 <= model.maxTourLength2 ) 
														&& ( r2 > 2 || route2.cost + changedDist_route2 <= model.maxTourLength2 ) )
												{
													//safe this move as a candidate move
													CROSS_option opt = new  CROSS_option( saving_CROSS, r1, r2, set1_start, set1_end, set2_start, set2end, -1);
													opt.capacityChangesRoute1 = set2_cap - set1_cap;
													opt.capacityChangesRoute2 = set1_cap - set2_cap;
													opt.costChangesRoute1 = changedDist_route1;
													opt.costChangesRoute2 = changedDist_route2;
													opt.directionSet1 = set1_direction;
													opt.directionSet2 = set2_direction;
													cross_options.add( opt );
												}//saving
											}//feasible
											//extend set2
											set2end += set2_direction;
											set2_nodes.add( route2.nodes.get(set2end).index );
											set2_cap += route2.nodes.get(set2end).demand;
											for (int pr=0; pr< model.products; pr++)
												invNeeded2[pr] += route2.nodes.get(set2end).demandPerProduct[pr];
											set2_dist+= model.getDistance(route2.nodes.get(set2end).index, route2.nodes.get(set2end-set2_direction).index);
										}//extension set2
										//extend set1
										set1_end += set1_direction;
										set1_nodes.add( route1.nodes.get(set1_end).index );
										set1_cap += route1.nodes.get(set1_end).demand;
										for (int pr=0; pr < model.products; pr++)
											invNeeded1[pr] += route1.nodes.get(set1_end).demandPerProduct[pr];
										set1_dist+= model.getDistance(route1.nodes.get(set1_end).index, route1.nodes.get(set1_end-set1_direction).index);
									}//extension set1
								}//CROSS=exchange CASE
								
								//-----------------------------------------------------------------------------
								//Crossover move
								//-----------------------------------------------------------------------------
								if ( model.depots == 1)
								{
									int set1_cap		=0;
									int set2_cap		=0;
									int set1_start		=node1_neighbour.routePosition ;
									int set2_start		=node2_pos;
									int set1_end		=route1.length;
									int set1_dir		=1;
									if (edgeNode == 1)
									{
										set1_end		=1;
										set1_dir		=-1;
									}
									int set2end			=route2.length;
									int set2_dir		=1;
									if (possibleStarts.get(0).directionRoute2 == 1)
									{
										set2end			=1;
										set2_dir		=-1;
									}
									
									//cannot start from a depot node
									if ( set1_start != 0 && set2_start != 0)
									{
										//check capacity constraint
										//TODO same with distance
										for (int i=set1_start; i!=set1_end+set1_dir; i+=set1_dir)
											set1_cap	+= route1.nodes.get(i).demand;
										for (int i=set2_start; i!=set2end+set2_dir; i+=set2_dir)
											set2_cap	+= route2.nodes.get(i).demand;
											
										if ( route1.load - set1_cap + set2_cap <= model.capacityLimitVehicle && route2.load - set2_cap + set1_cap <= model.capacityLimitVehicle)
										{model.ev2++;
											//create the move
											CROSS_option opt = new  CROSS_option( possibleStarts.get(0).savings, r1, r2, set1_start, set1_end, set2_start, set2end, -1);
											opt.capacityChangesRoute1 = set2_cap - set1_cap;
											opt.capacityChangesRoute2 = set1_cap - set2_cap;
											opt.costChangesRoute1 = -possibleStarts.get(0).savings;
											opt.costChangesRoute2 = 0; //TODO this is wrong
											opt.directionSet1 = set1_dir;
											opt.directionSet2 = set2_dir;
											cross_options.add( opt );
										}
									}
								}//Crossover
								
								}//Cross moves
								
								//------------------------------------------------------------------------------------------
								//OR-opt move (one set is emty)
								//------------------------------------------------------------------------------------------
								else
								{
									double saving_start =  possibleStarts.get(0).savings;
									int set2_direction  =  possibleStarts.get(0).directionRoute2;
									int set1_direction = 1;
									if (edgeNode == 0)
										set1_direction = -1;
									if ( possibleStarts.get(0).setEmpty == 1 )
										set1_direction *= -1;
									else
										set2_direction *= -1;
									
									//case 1: nodes come from route1 into route2
									int set_end 			= node1_pos;
									int set_direction 		= set1_direction;
									Route route_receive 	= route2;
									Route route_give		= route1;
									int route_give_index    = r1;//MTVRP
									int route_receive_index = r2;
									int set_corner 			= node1_neighbour.index;
									int set_reconnect		= route2.nodes.get(node2_pos + set2_direction).index;
									//case 2: nodes come from route2 into route1
									if ( possibleStarts.get(0).setEmpty == 1 )
									{
										set_end 		= node2_pos;
										set_direction 	= set2_direction;
										route_receive 	= route1;
										route_give		= route2;
										route_give_index    = r2;
									    route_receive_index = r1;
										set_corner 		= route2.nodes.get(node2_pos - set2_direction).index;
										set_reconnect	= node1_neighbour.index;
									}
									int set_start = set_end;
									int set_length = 1;
								
									ArrayList<Integer> set_nodes = new ArrayList<Integer>();
									set_nodes.add( route_give.nodes.get(set_end).index );
								
									int set_cap = route_give.nodes.get(set_end).demand;
									int[] invNeeded = new int[model.products];
									for (int pr=0; pr < model.products; pr++)
										invNeeded[pr] = route_give.nodes.get(set_end).demandPerProduct[pr];																
									double dist_set = 0;
									
									//extend set until the depot node is reached, or constraints are violated								
									while (set_end != 0 && set_end != route_give.length+1 
									&& route_receive.load + set_cap <= model.capacityLimitVehicle)
									{
										//feasible in terms of capacity and inventory?
										if (model.inventoryOff || depot1==depot2 || model.hasSufficientInventory(route_receive.depot, invNeeded)) 
										if (set_length > 0)
										{
											//get all relevant distances
											double dist_set_end 		= model.getDistance(route_give.nodes.get(set_end).index, set_reconnect);
											double dist_closure 		= model.getDistance(route_give.nodes.get(set_end+set_direction).index, set_corner);
											double dist_set_end_break_p = model.getDistance2(route_give.nodes.get(set_end).index, route_give.nodes.get(set_end+set_direction).index);
											double dist_set_end_break 	= model.getDistance(route_give.nodes.get(set_end).index, route_give.nodes.get(set_end+set_direction).index);
											model.ev2++;
											//compute the changes to both routes in terms of length
											double changedDist_receive = 0;
											double changedDist_give = 0;
											
											if ( possibleStarts.get(0).setEmpty == 1 ) //nodes go into route1
											{
												changedDist_receive = dist_n1_n2 + dist_set + dist_set_end - dist_n1_break;
												changedDist_give	= dist_closure - possibleStarts.get(0).dist_n2_break_T - dist_set_end_break - dist_set;
											}
											else //nodes go into route2
											{
												changedDist_receive = dist_n1_n2 + dist_set + dist_set_end - possibleStarts.get(0).dist_n2_break_T;
												changedDist_give	= dist_closure - dist_n1_break  - dist_set_end_break - dist_set;
											}

											//compute the saving of the move
											double saving_end 		= dist_set_end_break_p - dist_closure - dist_set_end;
											double saving_CROSS		= saving_start + saving_end;
											
											//LRP
											 if (route_give.length == set_length)
                                                 saving_CROSS += model.cost_per_Route;
											
											//does the move result in a saving and is it feasible in terms of length
											if ( saving_CROSS >= 0
													&& route_give.cost + changedDist_give <= model.maxTourLength 
													&& route_receive.cost + changedDist_receive <= model.maxTourLength 
													&& ( route_give_index> 2 || route_give.cost + changedDist_give <= model.maxTourLength2 ) 
													&& ( route_receive_index > 2 || route_receive.cost + changedDist_receive <= model.maxTourLength2 ) )
											{
												CROSS_option opt = null;
												if ( possibleStarts.get(0).setEmpty == 1 )//nodes go into route1
												{
													//create the candidate move
													opt = new  CROSS_option( saving_CROSS, r1, r2, Math.max(node1_pos, node1_pos + set1_direction), 999999, set_start, set_end, Math.max(node1_pos, node1_pos + set1_direction));
													opt.capacityChangesRoute1 	= set_cap;
													opt.capacityChangesRoute2 	= -set_cap;
													opt.costChangesRoute1 		= changedDist_receive;
													opt.costChangesRoute2 		= changedDist_give;
												}
												else //into Route 2
												{
													//create the candidate move
													opt = new  CROSS_option( saving_CROSS, r1, r2, set_start, set_end, Math.max(node2_pos, node2_pos + set2_direction), 999999, Math.max(node2_pos, node2_pos + set2_direction));
													opt.capacityChangesRoute1 	= -set_cap;
													opt.capacityChangesRoute2 	= set_cap;
													opt.costChangesRoute1 		= changedDist_give;
													opt.costChangesRoute2 		= changedDist_receive;
												}
												//get the nodes that surround set1 and set2
												opt.directionSet1 = set1_direction;
												opt.directionSet2 = set2_direction;
												cross_options.add( opt );						
											}//saving
										}//feasible
										//extend the set
										set_length++;
										set_end += set_direction;
										set_nodes.add( route_give.nodes.get(set_end).index );
										set_cap += route_give.nodes.get(set_end).demand;
										for (int pr=0; pr< model.products; pr++)
											invNeeded[pr] += route_give.nodes.get(set_end).demandPerProduct[pr];
										dist_set+= model.getDistance(route_give.nodes.get(set_end).index, route_give.nodes.get(set_end-set_direction).index);//route_give.nodes.get(set_end).distances[route_give.nodes.get(set_end-set_direction).index];
									}//extend set
								}//else
							
							possibleStarts.remove(0);					
						}//for all possible starts
								
					}//for all nearest nodes
				}//if not depot
			}//start from both incident nodes of the considered edge
		
		//sort all CROSS options according to their savings
		Collections.sort( cross_options );
		
		//----------------------------------------------------
		//execute the candidate move with the most savings
		//----------------------------------------------------
		while  ( !cross_options.isEmpty() )
		{	
			//execute the CROSS operation with the most savings
			CROSS_option cross1 = cross_options.get( 0 );
			executeCROSS( cross1 );
			
			ArrayList<Integer> forbiddenRoutes = new ArrayList<Integer>();	
			forbiddenRoutes.add( cross1.route1 ); forbiddenRoutes.add( cross1.route2 );
			
			//remove all candidate moves with the same routes
			Iterator<CROSS_option> it2 = cross_options.iterator();
			while (it2.hasNext())
			{
				//TODO inventory check!
				CROSS_option nextElement = it2.next();
				int r1 = nextElement.route1;
				int r2 = nextElement.route2;
				if ( forbiddenRoutes.contains(r1) || forbiddenRoutes.contains(r2) )
					it2.remove();		
			}//while
		}//while
	}//function
	
	
	private  void executeCROSS( CROSS_option cross )
	{
		double before = model.computeSolutionCosts();
		model.ex2++;
		
		//-------------------------------------------------------------------------------------------
		//nodes at positions [set1_start, set1_end] are inserted in route 2
		//nodes at positions [set2_start, set2_end] are inserted in route 1
		//-------------------------------------------------------------------------------------------
		int set1_start 	= Math.min( cross.set1_start, cross.set1_end );
		int set1_end 	= Math.max( cross.set1_start, cross.set1_end );
		int set2_start 	= Math.min( cross.set2_start, cross.set2_end );
		int set2_end	= Math.max( cross.set2_start, cross.set2_end );
		Route route1 	= model.routes.get( cross.route1 );
		Route route2 	= model.routes.get( cross.route2 );		
		
		boolean reverseSets = false;
		if ( cross.directionSet1 != cross.directionSet2 ) 
			reverseSets = true;
		
		ArrayList<Node> set1_nodes = new ArrayList<Node>();
		ArrayList<Node> set2_nodes = new ArrayList<Node>();
		
		//1) remove nodes from both routes
		if ( set1_end < 99999 )//if the set is not empty
		for (int i=set1_start; i <= set1_end; i++)
		{
			//adapt inventory
			if (!model.inventoryOff)
			for (int p=0; p<model.products; p++)
			{
				model.remainingInventory[route1.depot][p] += route1.nodes.get(set1_start).demandPerProduct[p];
				model.remainingInventory[route2.depot][p] -= route1.nodes.get(set1_start).demandPerProduct[p];
			}
			//remove node
			set1_nodes.add( route1.nodes.get(set1_start) );
			route1.nodes.remove(set1_start);
			route1.length--;
		}
		if ( set2_end < 99999 )//if the set is not empty
		for (int i=set2_start; i <= set2_end; i++)
		{
			//adapt inventory
			if (!model.inventoryOff)
			for (int p=0; p<model.products; p++)
			{
				model.remainingInventory[route2.depot][p] += route2.nodes.get(set2_start).demandPerProduct[p];
				model.remainingInventory[route1.depot][p] -= route2.nodes.get(set2_start).demandPerProduct[p];
			}
			//remove node
			set2_nodes.add( route2.nodes.get(set2_start) );
			route2.nodes.remove(set2_start);
			route2.length--;
		}

		//1 b) reverse the sets of nodes if the direction in both routes differs
		if ( reverseSets )
		{
			Collections.reverse(set1_nodes);
			Collections.reverse(set2_nodes);
		}

		//2) reinsert nodes into the other route
		route1.nodes.addAll(set1_start, set2_nodes);
		route1.length += set2_nodes.size();
		route1.load += cross.capacityChangesRoute1;
		route1.cost += cross.costChangesRoute1;
		
		route2.nodes.addAll(set2_start, set1_nodes);
		route2.length += set1_nodes.size();
		route2.load += cross.capacityChangesRoute2;
		route2.cost += cross.costChangesRoute2;
		
		//System.out.println("OUT " + set1_nodes.size() + "-" + set2_nodes.size() + " CROSS" );
		
		//3) post-processing
		//a) prepare intra-route optimisation
		if (!model.routes_to_be_optimised.contains( cross.route1 ) )
			model.routes_to_be_optimised.add( cross.route1 );
		if (!model.routes_to_be_optimised.contains( cross.route2 ) )
			model.routes_to_be_optimised.add( cross.route2 );
		if (!model.routes_changed_during_perturbation.contains( cross.route1 ) )
			model.routes_changed_during_perturbation.add( cross.route1 );
		if (!model.routes_changed_during_perturbation.contains( cross.route2 ) )
			model.routes_changed_during_perturbation.add( cross.route2 );
		if (!model.routes_update_RC_variables.contains(cross.route1)) model.routes_update_RC_variables.add( cross.route1 );
		if (!model.routes_update_RC_variables.contains(cross.route2)) model.routes_update_RC_variables.add( cross.route2 );
		if (!model.routes_update_edge_values.contains(cross.route1)) model.routes_update_edge_values.add( cross.route1 );
		if (!model.routes_update_edge_values.contains(cross.route2)) model.routes_update_edge_values.add( cross.route2 );
		
		//b) update customers
		for (int i=1; i<route1.nodes.size()-1; i++)
		{
			route1.nodes.get(i).routePosition 	= i;
			route1.nodes.get(i).routeIndex 		= cross.route1;	
			route1.nodes.get(i).depot			= route1.depot;
		}
		for (int i=1; i<route2.nodes.size()-1; i++)
		{
			route2.nodes.get(i).routePosition 	= i;
			route2.nodes.get(i).routeIndex 		= cross.route2;		
			route2.nodes.get(i).depot			= route2.depot;
		}
		
		//testing
		//if (model.routes.get(saving.from).nodes.get(0).depot != model.routes.get(saving.from).nodes.get(model.routes.get(saving.from).length+1).depot)
			//System.out.println("ERROR - after CROSS one tour is associated with two depots");

		double after = model.computeSolutionCosts();
		double diff = after - before;
		//if (diff < cross.costChangesRoute1 + cross.costChangesRoute2 - 0.1 || diff > cross.costChangesRoute1 + cross.costChangesRoute2+0.1)
			//System.out.println("ERROR : CROSS operation did not result in expected savings");
		
		//model.checkSolution();

		if (!model.penalties_off && !model.heuristicStart)
			model.moves_CROSS += 1;		
	}
}
