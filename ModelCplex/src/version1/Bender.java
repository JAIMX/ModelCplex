package version1;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import ilog.concert.*;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.Param.Preprocessing;
import version1.Data.Edge;

public class Bender {

	/**
	 * Minimize Z=c'x+f'y s.t Ax+By>=b x>=0 y in Y
	 * 
	 * Notice that here we exchange variable x and y compared with class Data
	 */

	private double UB, LB;
	private double[] c, f, b;
	private double[][] A, B;
	private int numOfX, numOfY, numOfConstraint, numOfTruck,numOfCity,numOfedge;
	private double[] currentY;
	private double[] b_By;
	private double fy;

	private IloCplex BMP, subProblem, extremeRay;
	private IloIntVar[] Y1;
	private IloNumVar z1, z3, dummy;
	private IloRange[] range2;
	private IloRange changeConstraint3;
	private IloObjective obj2, obj1;

	private IloNumVar[] u2, u3;
	private double tolerance;
	private double[] optX, optY;
	private final double diff = Math.pow(10, -5);

	// for BMP
	private ArrayList<double[]> feasibleCut, optimalCut;
	// private int numOfFeasibleCut,numOfOptimalCut;
	private ArrayList<Edge> edgeSet;
	private final int T;
	private ArrayList<ArrayList<Integer>> distance;
	private ArrayList<ArrayList<Integer>> distanceReverse;
	private int[] truckStartNode;

	public Bender(Data data, double tolerance) {
		// TODO Auto-generated constructor stub
		c = data.getc();
		f = data.getf();
		b = data.getbb();
		numOfX = data.getNumOfy();
		numOfY = data.getNumOfx();
		A = data.getA();
		B = data.getB();
		numOfConstraint = A.length;
		currentY = data.generateInitialx();
		this.tolerance = tolerance;
		optY = new double[numOfY];
		numOfTruck = data.getNumOfTruck();
		this.edgeSet = data.getEdgeSet();
		T=data.getT();
		numOfCity=data.getNumOfCity();
		distance=data.getDistance();
		distanceReverse=data.getDistanceReverse();
		numOfedge=edgeSet.size();
		truckStartNode=data.getTruckStartNode();
	}

	public void initailizeModel() throws IloException {
		// BMP
		/**
		 * minimize z1 s.t. z1+(u'B-f')y>=u'b u'By>=u'b y in Y
		 */
		feasibleCut = new ArrayList<double[]>();
		optimalCut = new ArrayList<double[]>();

		// subProblem
		/**
		 * max (b-B*currentY)'u s.t. A'u<=c u>=0 we only need to exchange objective
		 */

		subProblem = new IloCplex();
		u2 = new IloNumVar[numOfConstraint];
		for (int i = 0; i < numOfConstraint; i++) {
			u2[i] = subProblem.numVar(0, Double.MAX_VALUE, "u" + i);
		}
		range2 = new IloRange[numOfX];

		for (int column = 0; column < numOfX; column++) {
			IloLinearNumExpr constraint = subProblem.linearNumExpr();
			for (int row = 0; row < numOfConstraint; row++) {
				constraint.addTerm(A[row][column], u2[row]);
			}
			range2[column] = subProblem.addLe(constraint, c[column]);
		}

		b_By = new double[numOfConstraint];
		for (int i = 0; i < numOfConstraint; i++) {
			double temp = 0;
			for (int j = 0; j < numOfY; j++) {
				temp += B[i][j] * currentY[j];
			}
			b_By[i] = b[i] - temp;
		}

		fy = 0;
		for (int i = 0; i < numOfY; i++) {
			fy += f[i] * currentY[i];
		}

		subProblem.setOut(null);
		subProblem.setParam(Preprocessing.Presolve, false);
		// subProblem.exportModel("subProblem.lp");
		// System.out.println("subProblem finished.");

		// extremeRay
		/**
		 * max dummy s.t.dummy=0 u'*b_By>diff A'u<=0 u>=0
		 * 
		 * we only need to add 2nd constraint when we generate extremeRay
		 */
		extremeRay = new IloCplex();
		u3 = new IloNumVar[numOfConstraint];
		for (int i = 0; i < numOfConstraint; i++) {
			u3[i] = extremeRay.numVar(0, Double.MAX_VALUE, "u" + i);
		}
		dummy = extremeRay.numVar(0, 0, "dummy");
		IloLinearNumExpr obj = extremeRay.linearNumExpr();
		obj.addTerm(1, dummy);
		extremeRay.addMaximize(obj);
		extremeRay.addEq(obj, 0);

		for (int column = 0; column < numOfX; column++) {
			IloLinearNumExpr constraint = extremeRay.linearNumExpr();

			for (int row = 0; row < numOfConstraint; row++) {
				constraint.addTerm(A[row][column], u3[row]);
			}
			extremeRay.addLe(constraint, 0);
		}
		extremeRay.setOut(null);
		extremeRay.setParam(Preprocessing.Presolve, false);
		// extremeRay.exportModel("extremeRay.lp");
		// System.out.println("extremeRay finished.");

	}

	public void opt() throws IloException, FileNotFoundException {
		// PrintWriter out=new PrintWriter("tempout.txt");

		UB = Double.MAX_VALUE;
		LB = Double.MIN_VALUE;
		int count = 0;

		while (UB - LB > tolerance && count != 1) {
			count++;
			System.out.print("Now the upper bound= " + UB + "   ");
			System.out.println("lower bound= " + LB);
			// System.out.print("currentY= ");
			// System.out.println(Arrays.toString(currentY));
			System.out.println();

			// out.print("Now the upper bound= " + UB+" ");
			// out.println("lower bound= " + LB);
			// out.print("currentY= ");
			// out.println(Arrays.toString(currentY));
			// out.println();

			IloLinearNumExpr obj = subProblem.linearNumExpr(fy);

			for (int i = 0; i < numOfConstraint; i++) {
				obj.addTerm(b_By[i], u2[i]);
			}
			obj2 = subProblem.addMaximize(obj);

			System.out.println("||----Now solve the subProblem----||");
			System.out.println();

			// out.println("||----Now solve the subProblem----||");
			// out.println();

			if (!subProblem.solve()) { // infeasible or unbounded
				if (subProblem.getStatus().toString() == "Infeasible") {
					System.out.println("-->The subProblem is infeasible");
					System.out.println("The problem is unbounded!");

					// out.println("-->The subProblem is infeasible");
					// out.println("The problem is unbounded!");
					return;
				} else {// unbounded
					System.out.println("-->The subProblem is unbounded, so add a feasible cut to BMP");
					// out.println("-->The subProblem is unbounded, so add a feasible cut to BMP");

					// add the changing constraint to model extremeRay
					IloLinearNumExpr changingConstraint = extremeRay.linearNumExpr();

					for (int i = 0; i < numOfConstraint; i++) {
						changingConstraint.addTerm(b_By[i], u3[i]);
					}

					/**
					 * There is still some problem about the generation of extreme ray
					 */
					// changeConstraint3=extremeRay.addEq(1, changingConstraint);
					changeConstraint3 = extremeRay.addLe(diff, changingConstraint);

					extremeRay.solve();

					double[] u = extremeRay.getValues(u3);
					double[] tempc = new double[numOfY + 1];
					for (int i = 0; i < numOfY; i++) {
						for (int mul = 0; mul < numOfConstraint; mul++) {
							tempc[i] += u[mul] * B[mul][i];
						}
						tempc[i] = tempc[i] - f[i];
					}

					for (int mul = 0; mul < numOfConstraint; mul++) {
						tempc[numOfY] += u[mul] * b[mul];
					}

					feasibleCut.add(tempc);

					System.out.println("The extremeRay problem  is " + extremeRay.getStatus());
					System.out.println("Add a feasible cut to BMP.");

					// out.println("The extremeRay problem is "+ extremeRay.getStatus());
					// out.println("Add a feasible cut to BMP.");
					extremeRay.remove(changeConstraint3);

				}
			} else {
				System.out.println("-->The subProblem has a optimal solution");
				// out.println("-->The subProblem has a optimal solution");

				// subProblem has the optimal solution
				if (subProblem.getObjValue() < UB) {
					UB = subProblem.getObjValue();
					optX = subProblem.getDuals(range2);
					System.arraycopy(currentY, 0, optY, 0, currentY.length);
					System.out.println("New better feasible solution generates:");
					// System.out.println("optX= " + Arrays.toString(optX));
					// System.out.println("optY= " + Arrays.toString(optY));
					System.out.println("optmalObj= " + UB);
					System.out.println();

					// out.println("New better feasible solution generates:");
					// out.println("optX= " + Arrays.toString(optX));
					// out.println("optY= " + Arrays.toString(optY));
					// out.println("optmalObj= " + UB);
					// out.println();
				}

				System.out.println("Add a optimal cut to BMP");
				// out.println("Add a optimal cut to BMP");

				double[] u = subProblem.getValues(u2);
				double[] tempc = new double[numOfY + 1];

				for (int i = 0; i < numOfY; i++) {
					for (int mul = 0; mul < numOfConstraint; mul++) {
						tempc[i] += u[mul] * B[mul][i];
					}
				}

				for (int mul = 0; mul < numOfConstraint; mul++) {
					tempc[numOfY] += u[mul] * b[mul];
				}

				optimalCut.add(tempc);

			}
			subProblem.remove(obj2);

			// solve BMP and produce new currentY, fy, b_By
			System.out.println("||----Now solve the BMP----||");
			BMP();
			// System.out.println(BMP.toString());
			// // test++;
			// // if(test==1) {
			// // BMP.exportModel("aa.lp");
			// // }
			// if (!BMP.solve()) {// infeasible or unbounded
			// if (BMP.getStatus().toString() == "Infeasible") {
			// System.out.println("-->The BMP is infeasible");
			// System.out.println("The problem is infeasible!");
			// return;
			// } else { // unbounded
			// /**
			// * This situation may not appear because the initial BMP generate a constraint
			// * for z, it always has an optimal solution or infeasible solution.
			// */
			// System.out.print("-->The BMP is unbounded,");
			// // currentY
			// BMP.remove(obj1);
			// IloObjective temp = BMP.addMaximize();
			// BMP.solve();
			// currentY = BMP.getValues(Y1);
			// System.out.println(" and the new currentY is" + Arrays.toString(currentY));
			//
			// // b_by
			// for (int i = 0; i < numOfConstraint; i++) {
			// double temp2 = 0;
			// for (int j = 0; j < numOfY; j++) {
			// temp2 += B[i][j] * currentY[j];
			// }
			// b_by[i] = b[i] - temp2;
			// }
			//
			// // fy
			// fy = 0;
			// for (int i = 0; i < numOfY; i++) {
			// fy += f[i] * currentY[i];
			// }
			//
			// BMP.remove(temp);
			// BMP.add(obj1);
			// }
			// } else {
			// System.out.print("-->The BMP is has an optimal solution, ");
			// currentY = BMP.getValues(Y1);
			// System.out.println("and the new currentY is" + Arrays.toString(currentY));
			//
			// // b_by
			// for (int i = 0; i < numOfConstraint; i++) {
			// double temp2 = 0;
			// for (int j = 0; j < numOfY; j++) {
			// temp2 += B[i][j] * currentY[j];
			// }
			// b_by[i] = b[i] - temp2;
			// }
			//
			// // fy
			// fy = 0;
			// for (int i = 0; i < numOfY; i++) {
			// fy += f[i] * currentY[i];
			// }
			//
			// LB = BMP.getObjValue();
			// System.out.println("The new LB is " + LB);
			// }

		}

		// out.close();

		// System.out.println();
		// System.out.println("The optimal objective = " + UB);
		// System.out.println("The optX= " + Arrays.toString(optX));
		// System.out.println("The optY= " + Arrays.toString(optY));
		// System.out.println("Now the LB= " + LB);

		// System.out.println(feasibleCut.size());
		// System.out.println(optimalCut.size());

	}

	public void BMP() throws IloException {
		// Initialize BMP model
		IloCplex BMP = new IloCplex();
		IloObjective obj = BMP.addMinimize();
		IloRange[] optConstraint = new IloRange[optimalCut.size()];
		IloRange[] feaConstraint = new IloRange[feasibleCut.size()];
		IloRange[] truckConstraint = new IloRange[numOfTruck];
		// artificial variable
		IloNumVar[] z = new IloNumVar[numOfTruck + 2];
		ArrayList<Path> pathSet = new ArrayList<Path>();
		;

		// optimal cut
		for (int i = 0; i < optimalCut.size(); i++) {
			optConstraint[i] = BMP.addRange(optimalCut.get(i)[numOfY], Double.MAX_VALUE);
		}


		// feasible cut
		for (int i = 0; i < feasibleCut.size(); i++) {
			feaConstraint[i] = BMP.addRange(feasibleCut.get(i)[numOfY], Double.MAX_VALUE);
		}

		// sum of variable k equal to 1
		for (int k = 0; k < numOfTruck; k++) {
			truckConstraint[k] = BMP.addRange(1, 1);
		}

		// z0
		IloColumn tempColumn = BMP.column(obj, 1);
		for (int i = 0; i < optimalCut.size(); i++) {
			tempColumn = tempColumn.and(BMP.column(optConstraint[i], 1));
		}

		/**
		 * Here is a special  condition, when there is no optimal cut in BMP, we set z0=0 to get dual variables.
		 */
		if(optimalCut.size()==0) {
			IloRange tempConstraint = BMP.addRange(0, 0);
			tempColumn=tempColumn.and(BMP.column(tempConstraint,1));
		}	
		z[0] = BMP.numVar(tempColumn, Double.MIN_VALUE, Double.MAX_VALUE);
		

		// z1
		tempColumn = BMP.column(obj, Double.MAX_VALUE / 10000);
		for (int i = 0; i < feasibleCut.size(); i++) {
			tempColumn = tempColumn.and(BMP.column(feaConstraint[i], 1));
		}
		z[1] = BMP.numVar(tempColumn, 0, Double.MAX_VALUE);

		// z2-zk+1
		for (int index = 2; index <= numOfTruck + 1; index++) {
			tempColumn = BMP.column(obj, Double.MAX_VALUE / 10000);
			tempColumn = tempColumn.and(BMP.column(truckConstraint[index - 2], 1));
			z[index] = BMP.numVar(tempColumn, 0, Double.MAX_VALUE);
		}

		// start branch and price
		Stack<Node> stack = new Stack<Node>();
		// record the branch information
		ArrayList<HashSet<Integer>> notCover;
		int[][] cover=new int[numOfTruck][T+2];
//		Cover = new ArrayList<HashSet<Integer>>();
		notCover = new ArrayList<HashSet<Integer>>();

		for (int k = 0; k < numOfTruck; k++) {
			HashSet<Integer> temp = new HashSet<Integer>();
//			Cover.add(temp);
//			temp = new HashSet<Integer>();
			notCover.add(temp);
		}
		for(int k=0;k<numOfTruck;k++) {
			for(int t=0;t<T+2;t++) {
				cover[k][t]=-1;
			}
		}

		// root node
		Node root = new Node();
		root.extractCol = new HashSet<IloNumVar>();
		root.addCol = new HashSet<IloNumVar>();
		root.branchTruck = -1;
		root.ifAddCol = false;
		stack.add(root);

		// start DFS
		while (stack.size() > 0) {
			Node currentNode = stack.peek();

			if (currentNode.ifAddCol == false) { // we don't extract the node

				if (currentNode.branchTruck >= 0) {
					// delete some path in the current model firstly
					if (currentNode.ifCover == true) {
//						Cover.get(csurrentNode.branchTruck).add(currentNode.branchEdge);
						int branchEdgeIndex=currentNode.branchEdge;
						int startNode=edgeSet.get(branchEdgeIndex).start;

						if(startNode<numOfCity*(T+1)) {
							int time=startNode%(T+1);
							cover[currentNode.branchTruck][time]=branchEdgeIndex;
						}else { //branch edge start from the origin(cover[truck][T+1])
							cover[currentNode.branchTruck][T+1]=branchEdgeIndex;
						}
						
						
					} else {
						notCover.get(currentNode.branchTruck).add(currentNode.branchEdge);
					}

					// extract column
					if (currentNode.ifCover == true) { // Xak=1
						for (int i = 0; i < pathSet.size(); i++) {
							Path currentPath = pathSet.get(i);

							if (currentPath.ifInModel == true && currentPath.truckIndex == currentNode.branchTruck) {
								for (int edgeIndex : currentPath.edgeIndexSet) {
									if (edgeIndex != currentNode.branchEdge) {
										int start = edgeSet.get(currentNode.branchEdge).start;
										int end = edgeSet.get(currentNode.branchEdge).end;

										if (edgeSet.get(edgeIndex).start == start
												|| edgeSet.get(edgeIndex).end == end) {
											// extract currentPath in the model, and record on
											// currentNode.extractCol,change currentPath.ifinmodel
											currentNode.extractCol.add(currentPath.column);
											BMP.remove(currentPath.column);
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
										// extract currentPath in the model, and record on currentNode.extractCol,change currentPath.ifinmodel
										currentNode.extractCol.add(currentPath.column);
										BMP.remove(currentPath.column);
										currentPath.ifInModel = false;
										break;
									}
								}
							}
						}

					}
				}
				
				BMP.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Primal);
				
				
				
				// add new column(subproblem,according to cover and not cover;find a new class Path)
				int count=0;
				boolean stopReason=false;
				while(count<100) {
					BMP.solve();
					
					/// FIND AND ADD A NEW SHORTEST PATH///
					double[] optPrice=BMP.getDuals(optConstraint);
					double[] feaPrice=BMP.getDuals(feaConstraint);
					double[] truckPrice=BMP.getDuals(truckConstraint);
					
					// if all of k trucks don't have negative reduced cost, then check=false
					boolean check=false;
					for(int k=0;k<numOfTruck;k++) {
						
						double[] dpFunction=new double[numOfTruck*(T+1)];
						int[] pathRecord=new int[numOfTruck*(T+1)];
						for(int i=0;i<dpFunction.length;i++) {
							dpFunction[i]=Double.MAX_VALUE;
						}
						
						int startNode=numOfTruck*(T+1)+truckStartNode[k];
						//original node
						if(cover[k][T+1]>0) {
							int edgeIndex=cover[k][T+1];
							int pointTo=edgeSet.get(edgeIndex).end;
							pathRecord[pointTo]=edgeSet.get(edgeIndex).start;
							//calculate cost
							int columnIndex=numOfedge*k+edgeIndex;
							double cost=0;
							
							for(int i=0;i<optimalCut.size();i++) {
								cost+=optPrice[i]*optimalCut.get(i)[columnIndex];
							}
							for(int i=0;i<feasibleCut.size();i++) {
								cost+=feaPrice[i]*feasibleCut.get(i)[columnIndex];
							}
							cost=-cost;
							dpFunction[pointTo]=cost;
						}else { // no restriction
							for(int edgeIndex:distance.get(startNode)) {
								if(!notCover.get(k).contains(edgeIndex)) {
									//calculate cost
									int columnIndex=numOfedge*k+edgeIndex;
									double cost=0;
									
									for(int i=0;i<optimalCut.size();i++) {
										cost+=optPrice[i]*optimalCut.get(i)[columnIndex];
									}
									for(int i=0;i<feasibleCut.size();i++) {
										cost+=feaPrice[i]*feasibleCut.get(i)[columnIndex];
									}
									cost=-cost;
									
									int pointTo=edgeSet.get(edgeIndex).end;
									if(cost<dpFunction[pointTo]) {
										dpFunction[pointTo]=cost;
										pathRecord[pointTo]=startNode;
									}
								}
							}
						}
					}
					
					
				}
				
				
				
			}
		}

	}

	private class Path {
		IloNumVar column;
		HashSet<Integer> edgeIndexSet;
		int truckIndex;

		boolean ifInModel;
	}

	private class Node {
		
		HashSet<IloNumVar> extractCol, addCol;
		int branchTruck, branchEdge;
		boolean ifCover;
		
		boolean ifAddCol;
	}

	public static void main(String[] args) throws IOException, IloException {
		Data data = new Data();
		data.readData("out2.txt");
		data.graphTransfer();
		data.matrixGenerator();
		data.generateInitialx();

		Bender bender = new Bender(data, 0);
		bender.initailizeModel();
		bender.opt();

	}

}
