package version6;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import ilog.concert.*;
import ilog.cplex.*;
import version6.Data.Edge;

public class Bender {

	/**
	 * Minimize Z=c'*y+f'x s.t A*y+B*x>=b y>=0 X in X
	 */

	static final double FUZZ = 1.0e-7;

	IloCplex master;
	IloCplex sub;
	IloNumVar[] y;
	// IloIntVar[] x;
	// IloNumVar z;
	double[] c, f, b;
	double[][] A, B;
	int numX, numY, numConstraint;
	double[] xValues, yValues; // subproblem y values ||master x values
	IloRange[] subConstraint;
	// HashMap<IloConstraint, Integer> rhs; // here we use rhs to record the index
	// of constraints in sub
	Data data;
	double tolerance, UB, LB, zMaster;

	// for master
	ArrayList<double[]> feasibleCut, optimalCut;
	double[] optPrice;
	double[] feaPrice;
	double[] truckPrice;
	double[] initialX;

	int T, numOfCity, numOfedge, numOfTruck;
	private ArrayList<Edge> edgeSet;
	private ArrayList<ArrayList<Integer>> distance;

	public Bender(double[] c, double[] f, double[] b, double[][] A, double[][] B, Data data, double tolerance,
			double[] initialX) throws IloException {

		numY = c.length;
		numX = f.length;
		numConstraint = A.length;
		this.data = data;
		this.tolerance = tolerance;

		this.c = c;
		this.f = f;
		this.b = b;

		this.A = A;
		this.B = B;

		this.initialX = initialX;
		System.out.println("intialX is "+Arrays.toString(initialX));

		T = data.T;
		numOfCity = data.numberOfCities;
		numOfedge = data.edgeSet.size();
		this.edgeSet = data.edgeSet;
		this.distance = data.distance;
		numOfTruck = data.numberOfTrucks;

		// rhs = new HashMap<IloConstraint, Integer>();

		// set up the master problem(which initially has no constraint)

		BuildMaster();
		// master.setOut(null);

		// set up the subproblem
		sub = new IloCplex();
		BuildSub();

		solve();

	}

	public void BuildSub() throws IloException {

		y = new IloNumVar[numY];
		for (int p = 0; p < data.numberOfDemandPair; p++) {
			for (int edgeIndex = 0; edgeIndex < data.numOfEdge12; edgeIndex++) {
				Edge tempedge = data.edgeSet.get(edgeIndex);
				y[p * data.numOfEdge12 + edgeIndex] = sub.numVar(0, Double.MAX_VALUE,
						"y" + tempedge.u + "," + tempedge.t1 + "," + tempedge.v + "," + tempedge.t2 + "," + p);
			}
		}
		// y = sub.numVarArray(numY, 0, Double.MAX_VALUE);

		sub.addMinimize(sub.scalProd(c, y));
		// constrain to be satisfied -- record the constraints for use later
		subConstraint = new IloRange[numConstraint];
		IloLinearNumExpr expr = sub.linearNumExpr();

		// add constraints (initially all zero, which makes the subproblem
		// infeasible)
		// -- record the constraints for use later
		// -- also map each constraint to the corresponding binary variable in the
		// master (for decoding Farkas certificates)
		for (int row = 0; row < numConstraint; row++) {
			expr.clear();
			expr = sub.scalProd(A[row], y);
			subConstraint[row] = sub.addGe(expr, 0, "subConstraint_" + row);
			// rhs.put(subConstraint[row], master.diff(b[row], master.scalProd(B[row], x)));
			// rhs.put(subConstraint[row], row);
		}

		// System.out.println("initial subproblem is "+sub.toString());
		// disable presolving of the subproblem (if the presolver realizes the
		// subproblem is infeasible, we do not get a dual ray)
		sub.setParam(IloCplex.BooleanParam.PreInd, false);
		// force use of the dual simplex algorithm to get a Farkas certificate
		sub.setParam(IloCplex.IntParam.RootAlg, 2);
		// suppress subproblem output to reduce clutter
		sub.setOut(null);
	}

	public void BuildMaster() throws IloException {
		// x = new IloIntVar[numX];
		// for (int k = 0; k < data.numberOfTrucks; k++) {
		// for (int edgeIndex = 0; edgeIndex < data.edgeSet.size(); edgeIndex++) {
		// Edge tempedge = data.edgeSet.get(edgeIndex);
		// x[k * data.edgeSet.size() + edgeIndex] = master.intVar(0, 1,
		// "x" + tempedge.u + "," + tempedge.t1 + "," + tempedge.v + "," + tempedge.t2 +
		// "," + k);
		// }
		// }
		//
		// // x = master.intVarArray(numX, 0, 1);
		// z = master.numVar(Double.MIN_VALUE, Double.MAX_VALUE, "ArificialVar");
		// master.addMinimize(master.sum(z, master.scalProd(f, x)), "Obj");

		// record the parameter of constraints in master
		feasibleCut = new ArrayList<double[]>();
		optimalCut = new ArrayList<double[]>();

	}

	public final void solve() throws IloException {
		// master.solve();
		// zMaster = master.getValue(z);
		// xValues = master.getValues(x);

		xValues = initialX;
		

		double tempConst = 0;
		for (int i = 0; i < numX; i++) {
			tempConst += f[i] * xValues[i];
		}

		UB = Double.MAX_VALUE;
		LB = Double.MIN_VALUE;

		while (UB - LB > tolerance) {

			System.out.println("Now the upper bound= " + UB);
			System.out.println("lower bound= " + LB);
			// System.out.print("currentX= ");
			// System.out.println(Arrays.toString(xValues));
			System.out.println();

			// set the supply constraint right-hand sides in the subproblem
			for (int row = 0; row < numConstraint; row++) {
				double temp = 0;
				for (int i = 0; i < numX; i++) {
					temp += B[row][i] * xValues[i];
				}
				subConstraint[row].setLB(b[row] - temp);
			}

			// solve the subproblem
			// System.out.println("||----Now solve the subProblem----||");
			// System.out.println();

			sub.solve();
			IloCplex.Status status = sub.getStatus();
			// IloNumExpr expr = master.numExpr();

			// System.out.println("The subproblem is "+status.toString());

			if (status == IloCplex.Status.Infeasible) {
				// subproblem is infeasible -- add a feasibility cut
				// first step: get a Farkas certificate, corresponding to a dual ray
				// along which the dual is unbounded

				IloConstraint[] constraints = new IloConstraint[numConstraint];
				double[] coefficients = new double[numConstraint];
				sub.dualFarkas(constraints, coefficients);

				double[] mu = coefficients;
				double[] tempFeaCut = new double[numX + 1];

				// add a feasibility cut
				double temp = 0;
				for (int i = 0; i < numConstraint; i++) {
					temp += mu[i] * b[i];
				}
				tempFeaCut[numX] = temp;

				for (int i = 0; i < numX; i++) {
					double para = 0;
					for (int j = 0; j < numConstraint; j++) {
						para += mu[j] * B[j][i];
					}
					// para = -para;
					// para=para-f[i];

					tempFeaCut[i] = para;
					// expr = master.sum(expr, master.prod(para, x[i]));
				}

				feasibleCut.add(tempFeaCut);
				// IloConstraint r = master.addGe(0, expr);
				// System.out.println("\n>>> Adding feasibility cut: " + r + "\n");
				System.out.println("\n>>> Adding feasibility cut: " + "\n");
			} else if (status == IloCplex.Status.Optimal) {

				if (sub.getObjValue() + tempConst + FUZZ < UB) {
					UB = sub.getObjValue() + tempConst;
					yValues = sub.getValues(y);
					System.out.println("Find a  better solution, and updata UB= "+UB);
					System.out.println();
				}

				// add an optimality cut
				double[] mu = sub.getDuals(subConstraint);
				double[] tempOptCut = new double[numX + 1];

				double temp = 0;
				for (int i = 0; i < numConstraint; i++) {
					temp += mu[i] * b[i];
				}
				tempOptCut[numX] = temp;

				for (int i = 0; i < numX; i++) {

					double para = 0;
					for (int j = 0; j < numConstraint; j++) {
						para += mu[j] * B[j][i];
					}
					// para = -para;
					para = para - f[i];

					tempOptCut[i] = para;
					// expr = master.sum(expr, master.prod(para, x[i]));
				}

				optimalCut.add(tempOptCut);

				// IloConstraint r = master.addGe(z, expr);
				// System.out.println("\n>>> Adding optimality cut: " + r + "\n");
				System.out.println("\n>>> Adding optimality cut: " + "\n");

			} else {
				// unexpected status -- report but do nothing
				System.err.println("\n!!! Unexpected subproblem solution status: " + status + "\n");
			}

			// update xValues and LB
			BMP();
			// xValues=BMP();
			//// xValues = master.getValues(x);
			tempConst = 0;
			for (int i = 0; i < numX; i++) {
				tempConst += f[i] * xValues[i];
			}
			// // LB = master.getValue(z)+tempConst;
			// LB = master.getObjValue();

		}

		System.out.println("Now the upper bound= " + UB);
		System.out.println("lower bound= " + LB);
		System.out.println();

		// if (master.solve()) {
		// System.out.println("optimal obj= " + master.getObjValue());
		// // double[] xValues = master.getValues(x);
		// // System.out.println("x= " + Arrays.toString(xValues));
		// // System.out.println("y= " + Arrays.toString(yValues));
		// } else {
		// System.out.println("The last master's status is " +
		// master.getStatus().toString());
		// }

	}
	
//	public String getName(int index) {
//		String string = "";
//		// Vst
//		if (index < numOfCity * (T + 1)) {
//			int node = index / (T + 1);
//			int time = index % (T + 1);
//			string = "v" + node + "," + time;
//		} else {
//			if (index < numOfCity * (T + 2)) {
//				int node = index % (numOfCity * (T + 1));
//				string = "o" + node;
//			} else {
//				int node = index % (numOfCity * (T + 2));
//				string = "d" + node;
//			}
//		}
//
//		return string;
//	}
	
	//name of X variables
	public String getName(int index) {
		String string = "";
		// k
		int k=index/numOfedge;
		int edgeIndex=index%numOfedge;
		Edge edge=edgeSet.get(edgeIndex);
		
		string="X"+edge.start+","+edge.end+","+k;

		return string;
	}

	private class Path {
		IloNumVar column;
		HashSet<Integer> edgeIndexSet;
		int truckIndex;

		boolean ifInModel;
	}

	private class Node {

		HashSet<Path> extractCol, addCol;
		int branchTruck, branchEdge;
		boolean ifCover;

		boolean ifAddCol;
	}

	// update xValues and LB
	public void BMP() throws IloException {

		double masterUB = Double.MAX_VALUE;
		double masterLB = Double.MIN_VALUE;
		double[] optSolution = new double[numX];
		double[] tempSolution = new double[numX];
		optSolution[0] = -1;

		// Initialize BMP model
		master = new IloCplex();
		IloObjective obj = master.addMinimize();
		IloRange[] optConstraint = new IloRange[optimalCut.size()];
		IloRange[] feaConstraint = new IloRange[feasibleCut.size()];
		IloRange[] truckConstraint = new IloRange[data.numberOfTrucks];

		// artificial variable
		IloNumVar[] z = new IloNumVar[data.numberOfTrucks + 2];
		ArrayList<Path> pathSet = new ArrayList<Path>();

		// optimal cut
		for (int i = 0; i < optimalCut.size(); i++) {
			optConstraint[i] = master.addRange(optimalCut.get(i)[numX], Double.MAX_VALUE);
			System.out.println("#" + i + " optimalCut is " + Arrays.toString(optimalCut.get(i)));
			
			for(int j=0;j<numX;j++) {
				if(Math.abs(optimalCut.get(i)[j])>0) {
					System.out.print(optimalCut.get(i)[j]+" "+getName(j)+", ");
				}
			}
			System.out.println();
		}

		// feasible cut
		for (int i = 0; i < feasibleCut.size(); i++) {
			feaConstraint[i] = master.addRange(feasibleCut.get(i)[numX], Double.MAX_VALUE);
			System.out.println("#" + i + " feasibleCut is " + Arrays.toString(feasibleCut.get(i)));
			///----------------------------------------------------check repeat feasible cut problem-------------------------------///
			for(int j=0;j<numX;j++) {
				if(Math.abs(feasibleCut.get(i)[j])>0) {
					System.out.print(feasibleCut.get(i)[j]+" "+getName(j)+", ");
				}
			}
			System.out.println();
		}
		


		// sum of variable k equal to 1
		for (int k = 0; k < data.numberOfTrucks; k++) {
			truckConstraint[k] = master.addRange(1, 1);
		}

		// z0
		IloColumn tempColumn = master.column(obj, 1);
		for (int i = 0; i < optimalCut.size(); i++) {
			tempColumn = tempColumn.and(master.column(optConstraint[i], 1));
		}

		// /**
		// * Here is a special condition, when there is no optimal cut in BMP, we set
		// z0=0
		// * to get dual variables.
		// */
		// if (optimalCut.size() == 0) {
		// IloRange tempConstraint = BMP.addRange(0, 0);
		// tempColumn = tempColumn.and(BMP.column(tempConstraint, 1));
		// }

		z[0] = master.numVar(tempColumn, 0, Double.MAX_VALUE);

		// z1
		// tempColumn = BMP.column(obj, Double.MAX_VALUE / 10000);
		tempColumn = master.column(obj, 100000000);
		for (int i = 0; i < feasibleCut.size(); i++) {
			tempColumn = tempColumn.and(master.column(feaConstraint[i], 1));
		}
		z[1] = master.numVar(tempColumn, 0, Double.MAX_VALUE);

		// z2-zk+1
		for (int index = 2; index <= data.numberOfTrucks + 1; index++) {
			// tempColumn = BMP.column(obj, Double.MAX_VALUE / 10000);
			tempColumn = master.column(obj, 100000000);
			tempColumn = tempColumn.and(master.column(truckConstraint[index - 2], 1));
			z[index] = master.numVar(tempColumn, 0, Double.MAX_VALUE);
		}

		/// ---------------------------Start branch and
		/// price--------------------------------///
		Stack<Node> stack = new Stack<Node>();

		// record the branch information
		ArrayList<HashSet<Integer>> notCover;
		int[][] cover = new int[data.numberOfTrucks][data.T + 2];
		notCover = new ArrayList<HashSet<Integer>>();

		for (int k = 0; k < data.numberOfTrucks; k++) {
			HashSet<Integer> temp = new HashSet<Integer>();
			notCover.add(temp);
		}
		for (int k = 0; k < data.numberOfTrucks; k++) {
			for (int t = 0; t < data.T + 2; t++) {
				cover[k][t] = -1;
			}
		}

		// root node
		Node root = new Node();
		root.extractCol = new HashSet<Path>();
		root.addCol = new HashSet<Path>();
		root.branchTruck = -1;
		root.ifAddCol = false;
		stack.add(root);

		// start DFS
		while (stack.size() > 0) {
			System.out.println("Now stack.size= " + stack.size());
//			System.out.println("Now the master model is " + master.toString());
			System.out.println();

			Node currentNode = stack.peek();

			if (currentNode.ifAddCol == false) {// we don't extract the node

				System.out.println("We don't extract the current node");

				if (currentNode.branchTruck >= 0) {
					// delete some path in the current model firstly
					System.out.println("1-Deal with the branch of current node");
					if (currentNode.ifCover == true) {
						// Cover.get(csurrentNode.branchTruck).add(currentNode.branchEdge);
						int branchEdgeIndex = currentNode.branchEdge;
						int startNode = data.edgeSet.get(branchEdgeIndex).start;

						if (startNode < data.numberOfCities * (data.T + 1)) {
							int time = startNode % (data.T + 1);
							cover[currentNode.branchTruck][time] = branchEdgeIndex;
						} else { // branch edge start from the origin(cover[truck][T+1])
							cover[currentNode.branchTruck][data.T + 1] = branchEdgeIndex;
						}

					} else {
						notCover.get(currentNode.branchTruck).add(currentNode.branchEdge);
					}

					System.out.println("2-Extract some column in BMP acoording to the branch");
					// extract column
					if (currentNode.ifCover == true) { // Xak=1
						for (int i = 0; i < pathSet.size(); i++) {
							Path currentPath = pathSet.get(i);

							if (currentPath.ifInModel == true && currentPath.truckIndex == currentNode.branchTruck) {
								for (int edgeIndex : currentPath.edgeIndexSet) {
									if (edgeIndex != currentNode.branchEdge) {
										int start = data.edgeSet.get(currentNode.branchEdge).start;
										int end = data.edgeSet.get(currentNode.branchEdge).end;

										if (data.edgeSet.get(edgeIndex).start == start
												|| data.edgeSet.get(edgeIndex).end == end) {
											// extract currentPath in the model, and record on
											// currentNode.extractCol,change currentPath.ifinmodel
											currentNode.extractCol.add(currentPath);
											master.remove(currentPath.column);
											currentPath.ifInModel = false;
											break;
										}
									}
								}
							}
						}
					} else { // Xak=0

						for (int i = 0; i < pathSet.size(); i++) {
							Path currentPath = pathSet.get(i);

							if (currentPath.ifInModel == true && currentPath.truckIndex == currentNode.branchTruck) {
								for (int edgeIndex : currentPath.edgeIndexSet) {
									if (edgeIndex == currentNode.branchEdge) {
										// extract currentPath in the model, and record on currentNode.extractCol,change
										// currentPath.ifinmodel
										currentNode.extractCol.add(currentPath);
										master.remove(currentPath.column);
										currentPath.ifInModel = false;
										break;
									}
								}
							}
						}

					}
				}

				master.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Primal);

				System.out.println("3-Start to add some path to BMP");
				// add new column(subproblem,according to cover and not cover;find a new class
				// Path)
				int count = 0;
				boolean check = true;
				boolean ifAllZero = true;

				/// -------------------------------------------------------------------------------------------///

				master.solve();
				if (master.getObjValue() < masterUB) {

//					if (master.getObjValue() > masterLB) {
//						masterLB = master.getObjValue();
//						System.out.println("masterLB update to " + masterLB);
//					}

					//
					// PrintWriter out= new PrintWriter("tempout.txt");
					// out.println("feasibleCut= ");
					//
					// for(int e=0;e<edgeSet.size();e++) {
					// for(int truck=0;truck<numOfTruck;truck++) {
					// out.print(feasibleCut.get(0)[truck*numOfTruck]+e+" ");
					// }
					// out.println();
					// }
					// out.println("feasibleCut[bi]= "+feasibleCut.get(0)[numOfTruck*numOfedge]);
					// out.close();

					int truckStartIndex = 0;
					for (;;) {

						master.exportModel("tempMaster.lp");
						System.out.println("tempMaster's outcome:");
						System.out.println("objective= "+master.getObjValue());
						System.out.println("z= "+Arrays.toString(master.getValues(z)));
						
						for (int pathIndex = 0; pathIndex < pathSet.size(); pathIndex++) {
							Path path = pathSet.get(pathIndex);
							if (path.ifInModel == true ) {
								System.out.println("path #"+pathIndex+"= "+master.getValue(path.column));
							}
						}
						
						/// FIND AND ADD A NEW SHORTEST PATH///
						optPrice = master.getDuals(optConstraint);
						feaPrice = master.getDuals(feaConstraint);
						truckPrice = master.getDuals(truckConstraint);

						System.out.print("optPrice= ");
						System.out.println(Arrays.toString(optPrice));
						System.out.print("feaPrice= ");
						System.out.println(Arrays.toString(feaPrice));
						System.out.print("truckPrice= ");
						System.out.println(Arrays.toString(truckPrice));

						// if all of k trucks don't have negative reduced cost, then check=false
						check = false;
						double totalPathCost = Double.MAX_VALUE;
						int lastEdgeIndex = -1;

						for (int kk = 0; kk < numOfTruck; kk++) {

							int k = (truckStartIndex + kk) % numOfTruck;

							// double[] dpFunction = new double[numOfTruck * (T + 1)];
							// int[] pathRecord = new int[numOfTruck * (T + 1)];
							double[] dpFunction = new double[numOfCity * (T + 1)];
							int[] pathRecord = new int[numOfCity * (T + 1)];
							for (int i = 0; i < dpFunction.length; i++) {
								dpFunction[i] = Double.MAX_VALUE;
							}

							int startNode = numOfCity * (T + 1) + data.truckStartNode[k];
							// original node
							if (cover[k][T + 1] > 0) {
								int edgeIndex = cover[k][T + 1];
								int pointTo = edgeSet.get(edgeIndex).end;
								pathRecord[pointTo] = edgeIndex;

								dpFunction[pointTo] = calculateCost(edgeIndex, k);
							} else { // no restriction
								for (int edgeIndex : distance.get(startNode)) {
									if (!notCover.get(k).contains(edgeIndex)) {

										double cost = calculateCost(edgeIndex, k);
										int pointTo = edgeSet.get(edgeIndex).end;
										if (cost < dpFunction[pointTo]) {
											dpFunction[pointTo] = cost;
											pathRecord[pointTo] = edgeIndex;
										}
									}
								}
							}

							// update t in[0,T]
							// Initialize
							int nextCoverTime = -1;
							int nextStartPoint = -1;
							for (int time = nextCoverTime + 1; time <= T; time++) {
								if (cover[k][time] > 0) {
									nextCoverTime = time;
									nextStartPoint = edgeSet.get(cover[k][time]).start;
									break;
								}
							}

							if (nextStartPoint < 0) {
								nextCoverTime = T + 1;
							}

							int currentTime = 0;
							// update dpFunction
							while (currentTime < T) {

								if (currentTime < nextCoverTime) {

									for (int node = 0; node < numOfCity; node++) {

										int nodeIndex = node * (T + 1) + currentTime;
										if (dpFunction[nodeIndex] < Double.MAX_VALUE / 100000) {

											for (int edgeIndex : distance.get(nodeIndex)) {
												Edge edge = edgeSet.get(edgeIndex);
												if ((edge.t2 < nextCoverTime
														|| (edge.t2 == nextCoverTime && edge.end == nextStartPoint))
														&& !notCover.get(k).contains(edgeIndex)) {

													double cost = calculateCost(edgeIndex, k);

													if (dpFunction[edge.end] > dpFunction[nodeIndex] + cost) {
														dpFunction[edge.end] = dpFunction[nodeIndex] + cost;
														pathRecord[edge.end] = edgeIndex;
													}

												}
											}
										}
									}

									currentTime++;

								} else {// currentTime=nextCoverTime
									int edgeIndex = cover[k][nextCoverTime];
									Edge edge = edgeSet.get(edgeIndex);

									double cost = calculateCost(edgeIndex, k);

									dpFunction[edge.end] = dpFunction[edge.start] + cost;
									pathRecord[edge.end] = edgeIndex;

									nextStartPoint = -1;
									for (int time = nextCoverTime + 1; time <= T; time++) {
										if (cover[k][time] > 0) {
											nextCoverTime = time;
											nextStartPoint = edgeSet.get(cover[k][time]).start;
											break;
										}
									}
									if (nextStartPoint < 0) {
										nextCoverTime = T + 1;
									}

									currentTime = edge.t2;

								}

							}

							// destination
							if (cover[k][T] > 0) {
								int edgeIndex = cover[k][T];
								Edge edge = edgeSet.get(edgeIndex);
								double cost = calculateCost(edgeIndex, k);

								totalPathCost = dpFunction[edge.start] + cost;
								lastEdgeIndex = edgeIndex;

							} else {
								for (int node = 0; node < numOfCity; node++) {
									int nodeIndex = node * (T + 1) + T;
									if (dpFunction[nodeIndex] < Double.MAX_VALUE / 100000) {
										for (int edgeIndex : distance.get(nodeIndex)) {
											if (edgeSet.get(edgeIndex).end == numOfCity*(T+2)+data.truckStartNode[k]) {

												double cost = calculateCost(edgeIndex, k);
												if (totalPathCost > dpFunction[nodeIndex] + cost) {
													totalPathCost = dpFunction[nodeIndex] + cost;
													lastEdgeIndex = edgeIndex;
												}
											}
										}
									}
								}
							}

							double pik = master.getDual(truckConstraint[k]);
							if (totalPathCost - pik < -FUZZ) { // add this path as column of BMP

								// a new path
								Path newPath = new Path();
								newPath.truckIndex = k;
								HashSet<Integer> edgeIndexSet = new HashSet<Integer>();

								// add all the edges to edgeIndexSet
								int currentEdgeIndex = lastEdgeIndex;
								edgeIndexSet.add(currentEdgeIndex);
								int tempNode = edgeSet.get(currentEdgeIndex).start;

								while (edgeSet.get(currentEdgeIndex).setIndex != 3) {
									currentEdgeIndex = pathRecord[tempNode];
									edgeIndexSet.add(currentEdgeIndex);
									tempNode = edgeSet.get(currentEdgeIndex).start;
								}

								// build up the column
								IloColumn addColumn = master.column(obj, 0);

								double para = 0;
								for (int i = 0; i < optimalCut.size(); i++) {
									para = 0;
									for (int edgeIndex : edgeIndexSet) {
										para += optimalCut.get(i)[k * numOfedge + edgeIndex];
									}
									addColumn = addColumn.and(master.column(optConstraint[i], para));

								}

								for (int i = 0; i < feasibleCut.size(); i++) {
									para = 0;
									for (int edgeIndex : edgeIndexSet) {
										para += feasibleCut.get(i)[k * numOfedge + edgeIndex];
									}
									addColumn = addColumn.and(master.column(feaConstraint[i], para));

								}

								addColumn = addColumn.and(master.column(truckConstraint[k], 1));

//								newPath.column = master.numVar(addColumn, 0, Double.MAX_VALUE);
								newPath.column = master.numVar(addColumn, 0, 1);
								newPath.ifInModel = true;
								newPath.edgeIndexSet = edgeIndexSet;
								pathSet.add(newPath);
								check = true;
								count++;
								System.out.println(
										"We add #" + count + " path to this node, and the path belong to truck " + k);
								
								System.out.println("Now the master problem is "+master.toString());

								/// --------------------------------------------------check shortest path
								/// problem-----------------///
								// if(count==1) {
								for (int edgeIndex : edgeIndexSet) {
									Edge temp = edgeSet.get(edgeIndex);
									System.out.println("start= " + temp.start + " at time " + temp.t1 + "     end= "
											+ temp.end + " at time " + temp.t2);
								}
								// }

								/// ----------------------------check
								/// end--------------------------------------------------------///
								truckStartIndex = k + 1;
								break;

							} else { // totalPathCost-pik>0

							}

						}

						master.solve();
						// check if all artificial vars are 0
						ifAllZero = true;
						for (int i = 1; i < numOfTruck + 2; i++) {
							System.out.print(master.getValue(z[i]) + " ");
							// System.out.println("Now i= " + i + " " + BMP.getValue(z[i]));
							if (Math.abs(master.getValue(z[i])) > FUZZ) {

								ifAllZero = false;
								break;
							}
						}

						if (!check || (count > 1000 && ifAllZero)) {
							master.exportModel("master.lp");
							break;
						}

					}
				}

				if (!check) {
					System.out.println("There is no path with negative reduced cost!");
				} else {
					System.out.println("There is still path with negative reduced cost.");
				}
				System.out.println("We add " + count + " paths to BMP.");
				// BMP.exportModel("tempout.lp");

				// // check if all artificial vars are 0
				// boolean ifAllZero = true;
				// for (int i = 1; i < numOfTruck + 2; i++) {
				// System.out.println("Now i= " + i + " " + BMP.getValue(z[i]));
				// if (Math.abs(BMP.getValue(z[i])) > RC_EPS) {
				// ifAllZero = false;
				// break;
				// }
				// }

				if (ifAllZero) {
					System.out.println("All artifacial vars equals 0");
				} else {
					System.out.println("Not all artifacial vars equals 0");
				}

				tempSolution = new double[numX];
				if ((check) || (!check && ifAllZero)) {
					// check if all vars are integeral
					System.out.println();
					System.out.println("Check if all vars are integeral?");
					boolean ifIntegeral = true;
					
					for (int k = 0; k < numOfTruck; k++) {
						for (int edge = 0; edge < numOfedge; edge++) {
//							Edge tempEdge=edgeSet.get(edge);
//							System.out.println(tempEdge.start+"->"+tempEdge.end);
							double sum = 0;

							for (int pathIndex = 0; pathIndex < pathSet.size(); pathIndex++) {
								Path path = pathSet.get(pathIndex);
								if (path.ifInModel == true && path.truckIndex == k
										&& path.edgeIndexSet.contains(edge)) {
									sum += master.getValue(path.column);
								}
							}
//							System.out.print(sum + " ");
							

							if (Math.abs(sum) < FUZZ || Math.abs(sum - 1) < FUZZ) {
								tempSolution[numOfedge * k + edge] = sum;
								if(sum>0) {
									System.out.println(getName(numOfedge * k + edge)+"="+sum);
								}
							} else {

								System.out.println("The answer is no, so we add a new branch!");
								Node node1 = new Node();
								HashSet<Path> extractCol = new HashSet<Path>();
								HashSet<Path> addCol = new HashSet<Path>();
								node1.extractCol = extractCol;
								node1.addCol = addCol;
								node1.branchTruck = k;
								node1.branchEdge = edge;
								node1.ifCover = false;
								node1.ifAddCol = false;
								stack.add(node1);

								Node node2 = new Node();
								extractCol = new HashSet<Path>();
								addCol = new HashSet<Path>();
								node2.extractCol = extractCol;
								node2.addCol = addCol;
								node2.branchTruck = k;
								node2.branchEdge = edge;
								node2.ifCover = true;
								node2.ifAddCol = false;
								stack.add(node2);

								// currentNode.ifAddCol=true;
								ifIntegeral = false;
								break;
							}
						}
						if (!ifIntegeral)
							break;

					}

					if (ifIntegeral) {
						// currentNode.ifAddCol = true;
						System.out.println("We find a new feasible solution.");
						if (master.getObjValue() < masterUB) {
							System.out.println("And the masterUB is updated");
							masterUB = master.getObjValue();
							System.arraycopy(tempSolution, 0, optSolution, 0, numX);
						}

					}

				}
				currentNode.ifAddCol = true;

				// end: we don't extract the node
			} else { // currentNode.ifAddCol == true, we extract the node
				currentNode = stack.pop();
				System.out.println("We should extract current node.");

				if (currentNode.branchTruck > 0) {

					// deal with extractCol,addCol
					for (Path path : currentNode.addCol) {
						master.remove(path.column);
						path.ifInModel = false;

					}

					System.out.println(currentNode.extractCol.size());
					for (Path path : currentNode.extractCol) {
						master.add(path.column);
						path.ifInModel = true;
					}

					// deal with cover and notCover
					if (currentNode.ifCover == false) {
						notCover.get(currentNode.branchTruck).remove(currentNode.branchEdge);
					} else {
						if (edgeSet.get(currentNode.branchEdge).setIndex == 3) {
							cover[currentNode.branchTruck][T + 1] = -1;
						} else {
							cover[currentNode.branchTruck][currentNode.branchEdge % (T + 1)] = -1;
						}
					}
				}

			}

		}

//		xValues = optSolution;
		System.arraycopy(optSolution, 0, xValues, 0, numX);
		System.out.println("xValues = "+Arrays.toString(xValues));
		
		for(int i=0;i<numX;i++) {
			if(xValues[i]>0) {
				
				System.out.println(getName(i)+"="+xValues[i]);
			}
		}
		LB = master.getObjValue();

		// if (Math.abs(optSolution[0]+1)<FUZZ) {
		// return Double.MIN_VALUE/10;
		// }else {
		// System.arraycopy(optSolution, 0, currentY, 0, numOfY);
		// return master.getObjValue();
		// }

	}

	public double calculateCost(int edgeIndex, int k) {
		// calculate cost
		int columnIndex = data.edgeSet.size() * k + edgeIndex;
		double cost = 0;

		for (int i = 0; i < optimalCut.size(); i++) {
			cost += optPrice[i] * optimalCut.get(i)[columnIndex];
		}
		for (int i = 0; i < feasibleCut.size(); i++) {
			cost += feaPrice[i] * feasibleCut.get(i)[columnIndex];
		}
		cost = -cost;
		return cost;
	}

	public static void main(String[] args) throws IOException, IloException {
		Data data = new Data();
		// data.readData("./data/temp.txt");
		// System.out.println("Read data done!");
		data.readData("./data/out_small.txt");
		// data.readData("./data/data1.txt");
		// data.readData("./data/data2.txt");
		data.graphTransfer();
		// System.out.println("Graph transfer done!");
		data.matrixGenerator();
		// System.out.println("MatrixGenerator done!");
		double tolerance = 0;

		Bender test = new Bender(data.c, data.f, data.bb, data.A, data.B, data, tolerance, data.generateInitialx());
	}

}
