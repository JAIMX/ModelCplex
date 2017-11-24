import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.util.*;

import ilog.concert.*;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;

public class MulticommodityFlowModel {

	private class demandPair {
		private int s;
		private int t;
		private int demandQuantity;

	}

	private HashMap<String, Integer> cityIndex;
	private HashMap<String, Integer> truckIndex;
	private int numberOfCities;
	private int numberOfTrucks;
	private ArrayList<demandPair> demandPairs;
	private int numberOfDemandPair;
	private double xx[];
	private double yy[];
	private double[][] length;
	private double[] openTime;
	private double[] closeTime;
	private double[] arrivalTime;
	private double[] processingTime;
	private int M;
	private double fixedCost;
	private double transportationCost;
	private int legLimit;
	private double distanceLimit;
	private double averageSpeed;
	private double drivingTimePerDay;
	private int[] truckCapacity;
	private int[] truckStartNode;
	private int T;
	private int[][] b;
	static double e = Math.pow(10, -2);

	private boolean[][] connect;

	private IloCplex cplex;
	private IloIntVar[][][] X;
	private IloNumVar[][][] y;

	private class Edge {
		int pointTo;
		double length;
		int setIndex;
		// hat A=1
		// AT=2
		// AO=3
		// AD=4
		int u;
		int v;
		int t1, t2;
	}

	private class Edge2 {
		int pointFrom;
		double length;
		int setIndex;
	}

	private ArrayList<ArrayList<Edge>> distance;
	private ArrayList<ArrayList<Edge2>> distanceReverse;

	public MulticommodityFlowModel() throws IOException {

		cityIndex = new HashMap<String, Integer>();
		truckIndex = new HashMap<String, Integer>();
		demandPairs = new ArrayList<demandPair>();

	}

	public void readData(String filename) throws IOException {

		Scanner in = new Scanner(Paths.get(filename));

		String temp = new String();
		// ---read "cities"---//
		temp = in.nextLine();
		assert (temp.substring(0, 3) == "City") : "Wrong cities";
		int lIndex = temp.indexOf("'", 0);
		int rIndex = temp.indexOf("'", lIndex + 1);
		int citiesIndex = 0;

		while (lIndex >= 0) {
			String city = temp.substring(lIndex + 1, rIndex);
			cityIndex.put(city, citiesIndex);
			citiesIndex++;
			lIndex = temp.indexOf("'", rIndex + 1);
			rIndex = temp.indexOf("'", lIndex + 1);
		}
		numberOfCities = cityIndex.size();

		// ---read "demandQuantity"---//
		temp = in.nextLine();
		assert (temp.substring(0, 5) == "demand") : "Wrong demandQuantity";
		temp = in.nextLine();

		while (!temp.equals("}")) {
			demandPair pair = new demandPair();

			lIndex = temp.indexOf("'", 0);
			rIndex = temp.indexOf("'", lIndex + 1);
			String city = temp.substring(lIndex + 1, rIndex);
			pair.s = cityIndex.get(city);
			// System.out.println("s="+pair.s);

			lIndex = temp.indexOf("'", rIndex + 1);
			rIndex = temp.indexOf("'", lIndex + 1);
			city = temp.substring(lIndex + 1, rIndex);
			pair.t = cityIndex.get(city);
			// System.out.println("t="+pair.t);

			int a = temp.lastIndexOf(": ");
			int demand = Integer.parseInt(temp.substring(a + 2, temp.length()));
			pair.demandQuantity = demand;
			demandPairs.add(pair);
			// System.out.println("demand="+pair.demandQuantity);

			temp = in.nextLine();
		}

		numberOfDemandPair = demandPairs.size();

		temp = in.nextLine();
		assert (temp.substring(0, 4) == "coord") : "Wrong coordinate";
		// ---read "coordinate"---//
		xx = new double[numberOfCities];
		yy = new double[numberOfCities];
		temp = in.nextLine();

		while (!temp.equals("}")) {
			lIndex = temp.indexOf("'", 0);
			rIndex = temp.indexOf("'", lIndex + 1);
			int city = Integer.parseInt(temp.substring(lIndex + 1, rIndex));

			lIndex = temp.indexOf("(", rIndex + 1);
			rIndex = temp.indexOf(",", lIndex + 1);
			xx[city] = Double.valueOf(temp.substring(lIndex + 1, rIndex));
			yy[city] = Double.valueOf(temp.substring(rIndex + 1, temp.length() - 1));
			// System.out.println(city);
			// System.out.println("x="+x[city]);
			// System.out.println("y="+y[city]);

			temp = in.nextLine();
		}

		length = new double[numberOfCities][numberOfCities];

		// Calculate distance
		for (int i = 0; i < numberOfCities; i++) {
			for (int j = i; j < numberOfCities; j++) {
				if (i == j) {
					length[i][j] = 0;
				} else {
					length[i][j] = Math.sqrt(Math.pow((xx[i] - xx[j]), 2) + Math.pow((yy[i] - yy[j]), 2));
					length[j][i] = length[i][j];
				}
			}
		}

		temp = in.nextLine();
		assert (temp.substring(0, 3) == "open") : "Wrong openTime";
		openTime = new double[numberOfCities];

		// ---read "openTime"---//
		lIndex = temp.indexOf("'", 0);
		rIndex = temp.indexOf("'", lIndex + 1);

		while (lIndex >= 0) {
			int city = cityIndex.get(temp.substring(lIndex + 1, rIndex));
			lIndex = temp.indexOf(":", rIndex + 1);
			rIndex = temp.indexOf(",", lIndex + 1);

			if (rIndex < 0) {
				openTime[city] = Double.valueOf(temp.substring(lIndex + 1, temp.length() - 1));
				break;
			} else {
				openTime[city] = Double.valueOf(temp.substring(lIndex + 1, rIndex));
			}

			lIndex = temp.indexOf("'", rIndex);
			rIndex = temp.indexOf("'", lIndex + 1);
		}

		temp = in.nextLine();
		assert (temp.substring(0, 4) == "close") : "Wrong closeTime";
		closeTime = new double[numberOfCities];

		// ---read "closeTime"---//
		lIndex = temp.indexOf("'", 0);
		rIndex = temp.indexOf("'", lIndex + 1);

		while (lIndex >= 0) {
			int city = cityIndex.get(temp.substring(lIndex + 1, rIndex));
			lIndex = temp.indexOf(":", rIndex + 1);
			rIndex = temp.indexOf(",", lIndex + 1);

			if (rIndex < 0) {
				closeTime[city] = Double.valueOf(temp.substring(lIndex + 1, temp.length() - 1));
				break;
			} else {
				closeTime[city] = Double.valueOf(temp.substring(lIndex + 1, rIndex));
			}

			lIndex = temp.indexOf("'", rIndex);
			rIndex = temp.indexOf("'", lIndex + 1);
		}

		temp = in.nextLine();
		assert (temp.substring(0, 6) == "arrival") : "Wrong arrivalTime";
		arrivalTime = new double[numberOfCities];

		// ---read "arrivalTime"---//
		lIndex = temp.indexOf("'", 0);
		rIndex = temp.indexOf("'", lIndex + 1);

		while (lIndex >= 0) {
			int city = cityIndex.get(temp.substring(lIndex + 1, rIndex));
			lIndex = temp.indexOf(":", rIndex + 1);
			rIndex = temp.indexOf(",", lIndex + 1);

			if (rIndex < 0) {
				arrivalTime[city] = Double.valueOf(temp.substring(lIndex + 1, temp.length() - 1));
				break;
			} else {
				arrivalTime[city] = Double.valueOf(temp.substring(lIndex + 1, rIndex));
			}

			lIndex = temp.indexOf("'", rIndex);
			rIndex = temp.indexOf("'", lIndex + 1);
		}
		T = (int) arrivalTime[0];

		temp = in.nextLine();
		assert (temp.substring(0, 6) == "process") : "Wrong processingTime";
		processingTime = new double[numberOfCities];

		// ---read "processingTime"---//
		lIndex = temp.indexOf("'", 0);
		rIndex = temp.indexOf("'", lIndex + 1);

		while (lIndex >= 0) {
			int city = cityIndex.get(temp.substring(lIndex + 1, rIndex));
			lIndex = temp.indexOf(":", rIndex + 1);
			rIndex = temp.indexOf(",", lIndex + 1);

			if (rIndex < 0) {
				processingTime[city] = Double.valueOf(temp.substring(lIndex + 1, temp.length() - 1));
				break;
			} else {
				processingTime[city] = Double.valueOf(temp.substring(lIndex + 1, rIndex));
			}

			lIndex = temp.indexOf("'", rIndex);
			rIndex = temp.indexOf("'", lIndex + 1);
		}

		// ---read single parameters---//
		// M
		temp = in.nextLine();
		lIndex = temp.indexOf("=", 0);
		M = Integer.parseInt(temp.substring(lIndex + 1, temp.length()));
		// fixedCost
		temp = in.nextLine();
		temp = in.nextLine();
		lIndex = temp.indexOf("=", 0);
		fixedCost = Double.valueOf(temp.substring(lIndex + 1, temp.length()));
		// transprotationCost
		temp = in.nextLine();
		lIndex = temp.indexOf("=", 0);
		transportationCost = Double.valueOf(temp.substring(lIndex + 1, temp.length()));
		// max number of legs per truck
		temp = in.nextLine();
		lIndex = temp.indexOf("=", 0);
		legLimit = Integer.parseInt(temp.substring(lIndex + 1, temp.length()));
		// max distance by a truck
		temp = in.nextLine();
		lIndex = temp.indexOf("=", 0);
		distanceLimit = Double.valueOf(temp.substring(lIndex + 1, temp.length()));
		// average speed
		temp = in.nextLine();
		lIndex = temp.indexOf("=", 0);
		averageSpeed = Double.valueOf(temp.substring(lIndex + 1, temp.length()));
		// driving time per day
		temp = in.nextLine();
		lIndex = temp.indexOf("=", 0);
		drivingTimePerDay = Double.valueOf(temp.substring(lIndex + 1, temp.length()));

		// ---read "trucks"---//
		temp = in.nextLine();
		assert (temp.substring(0, 5) == "trucks") : "Wrong trucks";
		temp = in.nextLine();
		int trucksIndex = 0;

		while (!temp.equals("}")) {
			lIndex = temp.indexOf("'", 0);
			rIndex = temp.indexOf("'", lIndex + 1);
			while (lIndex >= 0) {
				String truck = temp.substring(lIndex + 1, rIndex);
				truckIndex.put(truck, trucksIndex);
				// System.out.println(truck.toString());
				trucksIndex++;

				lIndex = temp.indexOf("'", rIndex + 1);
				rIndex = temp.indexOf("'", lIndex + 1);
			}
			temp = in.nextLine();
		}
		// System.out.println(truckIndex.size());
		numberOfTrucks = truckIndex.size();

		// ---read "truckCapacity"---//
		temp = in.nextLine();
		assert (temp.substring(0, 4) == "truck") : "Wrong trucksCapacity";
		truckCapacity = new int[numberOfTrucks];

		while (!temp.equals("}")) {
			lIndex = temp.indexOf("'", 0);
			rIndex = temp.indexOf("'", lIndex + 1);
			int index1 = temp.indexOf(":", rIndex + 1);
			int index2 = temp.indexOf(",", rIndex + 1);

			while (lIndex >= 0) {
				int capacity = Integer.valueOf(temp.substring(index1 + 1, index2));
				truckCapacity[truckIndex.get(temp.substring(lIndex + 1, rIndex))] = capacity;

				lIndex = temp.indexOf("'", index2 + 1);
				if (lIndex >= 0) {
					rIndex = temp.indexOf("'", lIndex + 1);
					index1 = temp.indexOf(":", rIndex + 1);
					index2 = temp.indexOf(",", rIndex + 1);
				}

			}
			temp = in.nextLine();
		}
		// ---read "truckStartNode"---//
		truckStartNode = new int[numberOfTrucks];

		temp = in.nextLine();
		assert (temp.substring(0, 4) == "truck") : "Wrong truckStartNode";
		temp = in.nextLine();

		while (!temp.equals("}")) {
			lIndex = temp.indexOf("'", 0);
			rIndex = temp.indexOf("'", lIndex + 1);
			int index1 = temp.indexOf(":", rIndex + 1);
			int index2 = temp.indexOf(",", rIndex + 1);

			while (lIndex >= 0) {

				int truckIndexTemp = truckIndex.get(temp.substring(lIndex + 1, rIndex));
				int cityIndexTemp = cityIndex.get(temp.substring(index1 + 1, index2));
				truckStartNode[truckIndexTemp] = cityIndexTemp;

				lIndex = temp.indexOf("'", index2 + 1);
				if (lIndex >= 0) {
					rIndex = temp.indexOf("'", lIndex + 1);
					index1 = temp.indexOf(":", rIndex + 1);
					index2 = temp.indexOf(",", rIndex + 1);
				}

			}
			temp = in.nextLine();
		}

	}

	/**
	 * Vst index:00,01,...,0T;n-1 0,n-1 1,...,n-1 T [0,n(T+1)-1] O
	 * index:n0,n1,...,n(n-1) [n(T+1),n(T+2)-1] D index:n+1 0,n+1 1,...,n+1(n-1)
	 * [n(T+2), n(T+3)-1]
	 */
	public void graphTransfer() {

		connect = new boolean[numberOfCities * (T + 3)][numberOfCities * (T + 3)];
		// calculate distance
		distance = new ArrayList<ArrayList<Edge>>();
		distanceReverse = new ArrayList<ArrayList<Edge2>>();

		// only record Vst and O
		for (int i = 0; i < numberOfCities * (T + 2); i++) {
			ArrayList<Edge> templist = new ArrayList<Edge>();
			distance.add(templist);

		}

		for (int i = 0; i < numberOfCities * (T + 3); i++) {
			ArrayList<Edge2> templist2 = new ArrayList<Edge2>();
			distanceReverse.add(templist2);
		}

		// add AT
		for (int node = 0; node < numberOfCities; node++) {
			for (int t = 0; t < T; t++) {
				Edge edge = new Edge();
				int nodeIndex = node * (T + 1) + t;
				edge.pointTo = nodeIndex + 1;
				edge.length = 0;
				edge.setIndex = 2;
				edge.u = node;
				edge.v = node;
				edge.t1 = t;
				edge.t2 = t + 1;
				distance.get(nodeIndex).add(edge);
				connect[nodeIndex][nodeIndex + 1] = true;

				Edge2 edge2 = new Edge2();
				edge2.pointFrom = nodeIndex;
				edge2.length = 0;
				edge2.setIndex = 2;
				distanceReverse.get(nodeIndex + 1).add(edge2);

			}
		}

		// add hat A
		for (int i = 0; i < numberOfCities; i++) {
			for (int j = 0; j < numberOfCities; j++) {
				if (i != j) {
					double time = length[i][j] / averageSpeed;
					int timeLength = (int) Math.ceil(time);

					int t = timeLength;
					int nodeIndex1 = i * (T + 1);
					int nodeIndex2 = j * (T + 1) + timeLength;

					while (t <= T) {
						Edge edge = new Edge();
						edge.pointTo = nodeIndex2;
						edge.length = length[i][j];
						edge.setIndex = 1;
						edge.u = i;
						edge.v = j;
						edge.t1 = t - timeLength;
						edge.t2 = t;
						distance.get(nodeIndex1).add(edge);
						connect[nodeIndex1][nodeIndex2] = true;

						Edge2 edge2 = new Edge2();
						edge2.pointFrom = nodeIndex1;
						edge2.length = edge.length;
						edge2.setIndex = 1;
						distanceReverse.get(edge.pointTo).add(edge2);

						t++;
						nodeIndex1++;
						nodeIndex2++;
					}
				}
			}
		}

		// add AO:(Ok,n0)
		for (int o = 0; o < numberOfCities; o++) {
			int oIndex = numberOfCities * (T + 1) + o;

			// for (int node = 0; node < numberOfCities; node++) {
			// int nodeIndex = node * (T + 1);
			// Edge edge = new Edge();
			// edge.pointTo = nodeIndex;
			// edge.length = length[o][node];
			// edge.setIndex = 3;
			// edge.u = o;
			// edge.v = node;
			// edge.t1 = -1;
			// edge.t2 = 0;
			// distance.get(oIndex).add(edge);
			// connect[oIndex][nodeIndex] = true;
			//
			// Edge2 edge2 = new Edge2();
			// edge2.pointFrom = oIndex;
			// edge2.length = edge.length;
			// edge2.setIndex = 3;
			// distanceReverse.get(edge.pointTo).add(edge2);
			//
			// }
			int nodeIndex = o * (T + 1);
			Edge edge = new Edge();
			edge.pointTo = nodeIndex;
			edge.length = 0;
			edge.setIndex = 3;
			edge.u = o;
			edge.v = o;
			edge.t1 = -1;
			edge.t2 = 0;
			distance.get(oIndex).add(edge);
			connect[oIndex][nodeIndex] = true;

			Edge2 edge2 = new Edge2();
			edge2.pointFrom = oIndex;
			edge2.length = 0;
			edge2.setIndex = 3;
			distanceReverse.get(edge.pointTo).add(edge2);
		}

		// add AD:(nT,Dk)
		for (int node = 0; node < numberOfCities; node++) {
			int nodeIndex = node * (T + 1) + T;
			// for (int d = 0; d < numberOfCities; d++) {
			// int dIndex = numberOfCities * (T + 2) + d;
			// Edge edge = new Edge();
			// edge.pointTo = dIndex;
			// edge.length = length[node][d];
			// edge.setIndex = 4;
			// edge.u = node;
			// edge.v = d;
			// edge.t1 = T;
			// edge.t2 = -1;
			// distance.get(nodeIndex).add(edge);
			// connect[nodeIndex][dIndex] = true;
			//
			// Edge2 edge2 = new Edge2();
			// edge2.pointFrom = nodeIndex;
			// edge2.length = edge.length;
			// edge2.setIndex = 4;
			// distanceReverse.get(edge.pointTo).add(edge2);
			// }

			int dIndex = numberOfCities * (T + 2) + node;
			Edge edge = new Edge();
			edge.pointTo = dIndex;
			edge.length = 0;
			edge.setIndex = 4;
			edge.u = node;
			edge.v = node;
			edge.t1 = T;
			edge.t2 = -2;
			distance.get(nodeIndex).add(edge);
			connect[nodeIndex][dIndex] = true;

			Edge2 edge2 = new Edge2();
			edge2.pointFrom = nodeIndex;
			edge2.length = 0;
			edge2.setIndex = 4;
			distanceReverse.get(edge.pointTo).add(edge2);
		}

		// set b
		b = new int[numberOfCities * (T + 1)][numberOfDemandPair];
		for (int p = 0; p < numberOfDemandPair; p++) {
			demandPair pair = demandPairs.get(p);
			b[pair.s * (T + 1)][p] = pair.demandQuantity;
			b[pair.t * (T + 1) + T][p] = -pair.demandQuantity;
		}

		// for(int i=0;i<b.length;i++) {
		// System.out.println(Arrays.toString(b[i]));
		// }

		// for(int i=0;i<numberOfCities;i++) {
		// System.out.println(Arrays.toString(length[i]));
		// }
		// for(int index=0;index<distance.size();index++) {
		// for(Edge e:distance.get(index)) {
		// System.out.println(e.u+","+e.t1+"->"+e.v+","+e.t2+": "+e.length+"
		// "+e.setIndex);
		// }
		// System.out.println();
		// }
	}

	public String getName(int index) {
		String string = "";
		// Vst
		if (index < numberOfCities * (T + 1)) {
			int node = index / (T + 1);
			int time = index % (T + 1);
			string = "v" + node + "," + time;
		} else {
			if (index < numberOfCities * (T + 2)) {
				int node = index % (numberOfCities * (T + 1));
				string = "o" + node;
			} else {
				int node = index % (numberOfCities * (T + 2));
				string = "d" + node;
			}
		}

		return string;
	}

	public void ModelBuilding() {
		try {
			cplex = new IloCplex();

			// ---decision variables---//
			int totalNumNode = numberOfCities * (T + 3);
			int numOfVst = numberOfCities * (T + 1);
			X = new IloIntVar[totalNumNode][totalNumNode][numberOfTrucks];
			y = new IloNumVar[numOfVst][numOfVst][numberOfDemandPair];

			for (int i = 0; i < totalNumNode; i++) {
				for (int j = 0; j < totalNumNode; j++) {
					for (int k = 0; k < numberOfTrucks; k++) {
						if (i != j) {
							X[i][j][k] = cplex.intVar(0, 1, "X" + i + "," + j + "," + k);
						} else {
							X[i][j][k] = cplex.intVar(0, 0, "X" + i + "," + j + "," + k);
						}
					}
				}
			}

			for (int i = 0; i < numOfVst; i++) {
				for (int j = 0; j < numOfVst; j++) {
					for (int p = 0; p < numberOfDemandPair; p++) {

						int s = demandPairs.get(p).s;
						int t = demandPairs.get(p).t;

						if (i != j) {
							y[i][j][p] = cplex.numVar(0, Double.MAX_VALUE, "y" + i + "," + j + "," + s + "," + t);
						} else {
							y[i][j][p] = cplex.numVar(0, 0, "y" + i + "," + j + "," + s + "," + t);
						}

					}
				}
			}

			// ---Objective---//
			IloLinearNumExpr obj = cplex.linearNumExpr();

			double parameter = fixedCost / (averageSpeed * drivingTimePerDay);

			for (int i = 0; i < numberOfCities * (T + 2); i++) {
				for (Edge e : distance.get(i)) {
					for (int k = 0; k < numberOfTrucks; k++) {
						obj.addTerm(parameter * e.length, X[i][e.pointTo][k]);
					}

				}
			}

			for (int i = 0; i < numberOfCities * (T + 1); i++) {
				for (Edge e : distance.get(i)) {
					if (e.pointTo < numberOfCities * (T + 1)) {
						for (int p = 0; p < numberOfDemandPair; p++) {
							obj.addTerm(transportationCost * e.length, y[i][e.pointTo][p]);
						}
					}

				}
			}

			cplex.addMinimize(obj);
			// System.out.println(cplex.getObjective().toString());

			// ---constraint 1-7---//
			for (int k = 0; k < numberOfTrucks; k++) {

				// ---constraint 1---//
				for (int node = 0; node < numberOfCities * (T + 1); node++) {

					IloLinearNumExpr constraint1 = cplex.linearNumExpr();
					for (Edge e : distance.get(node)) {
						constraint1.addTerm(1, X[node][e.pointTo][k]);
					}

					for (Edge2 e : distanceReverse.get(node)) {
						constraint1.addTerm(-1, X[e.pointFrom][node][k]);
					}

					// System.out.println(constraint1.toString());

					cplex.addEq(0, constraint1);
				}

				// ---constraint 2---//
				IloLinearNumExpr constraint2 = cplex.linearNumExpr();
				int ok = numberOfCities * (T + 1) + truckStartNode[k];
				for (Edge e : distance.get(ok)) {
					constraint2.addTerm(1, X[ok][e.pointTo][k]);
				}
				// System.out.println(constraint2.toString());
				cplex.addGe(1, constraint2);

				// ---constraint 3---//
				IloLinearNumExpr constraint3 = cplex.linearNumExpr();
				ok = numberOfCities * (T + 1) + truckStartNode[k];
				for (int o = numberOfCities * (T + 1); o < numberOfCities * (T + 2); o++) {
					if (o != ok) {
						for (Edge e : distance.get(o)) {
							constraint3.addTerm(1, X[o][e.pointTo][k]);
						}
					}
				}
				// System.out.println(constraint3.toString());
				cplex.addEq(0, constraint3);

				// // ---constraint 4---//
				// IloLinearNumExpr constraint4 = cplex.linearNumExpr();
				// int dk=numberOfCities*(T+2)+truckStartNode[k];
				// for(Edge2 e:distanceReverse.get(dk)) {
				// constraint4.addTerm(1, X[e.pointFrom][dk][k]);
				// }
				// cplex.addGe(1, constraint4);

				// ---constraint 5---//
				IloLinearNumExpr constraint5 = cplex.linearNumExpr();
				int dk = numberOfCities * (T + 2) + truckStartNode[k];
				for (int d = numberOfCities * (T + 2); d < numberOfCities * (T + 3); d++) {
					if (d != dk) {
						for (Edge2 e : distanceReverse.get(d)) {
							constraint5.addTerm(1, X[e.pointFrom][d][k]);
						}
					}
				}
				// System.out.println(constraint5);
				cplex.addEq(0, constraint5);

				// ---constraint 6---//
				IloLinearNumExpr constraint6 = cplex.linearNumExpr();
				for (int node = 0; node < distance.size(); node++) {
					for (Edge e : distance.get(node)) {
						if (e.setIndex == 1) {
							constraint6.addTerm(1, X[node][e.pointTo][k]);
						}

					}
				}

				// System.out.println(constraint6);
				cplex.addGe(legLimit, constraint6);

				// ---constraint 7---//
				IloLinearNumExpr constraint7 = cplex.linearNumExpr();
				for (int node = 0; node < distance.size(); node++) {
					for (Edge e : distance.get(node)) {
						constraint7.addTerm(e.length, X[node][e.pointTo][k]);
					}
				}

				// System.out.println(constraint7);
				cplex.addGe(distanceLimit, constraint7);

			}

			// ---constraint 8---//
			for (int p = 0; p < numberOfDemandPair; p++) {
				for (int i = 0; i < numberOfCities * (T + 1); i++) {
					IloLinearNumExpr constraint8 = cplex.linearNumExpr();
					for (Edge e : distance.get(i)) {
						if (e.pointTo < numberOfCities * (T + 2)) {
							constraint8.addTerm(1, y[i][e.pointTo][p]);
						}
					}

					for (Edge2 e : distanceReverse.get(i)) {
						if (e.pointFrom < numberOfCities * (T + 1)) {
							constraint8.addTerm(-1, y[e.pointFrom][i][p]);
						}
					}

					// System.out.println(constraint8);
					cplex.addEq(b[i][p], constraint8);
				}
			}

			// ---constraint 9---//
			for (int i = 0; i < distance.size(); i++) {
				for (Edge e : distance.get(i)) {
					if (e.setIndex == 1) {
						IloLinearNumExpr constraint9 = cplex.linearNumExpr();
						for (int p = 0; p < numberOfDemandPair; p++) {
							constraint9.addTerm(1, y[i][e.pointTo][p]);
						}

						for (int k = 0; k < numberOfTrucks; k++) {
							constraint9.addTerm(-truckCapacity[k], X[i][e.pointTo][k]);
						}

						// System.out.println(constraint9);
						cplex.addGe(0, constraint9);
					}
				}
			}

			// cplex.exportModel("MulticommodityFlowModel.lp");
			// cplex.setParam(IloCplex.Param.RootAlgorithm,
			// IloCplex.Algorithm.Primal);
			cplex.setParam(IloCplex.Param.Emphasis.Memory, true);
			cplex.setParam(IloCplex.IntParam.NodeFileInd, 2);
			cplex.setParam(IloCplex.IntParam.Threads, 1);
//			cplex.setParam(IloCplex.DoubleParam.EpGap, 0.025);
			// cplex.setParam(IloCplex.IntParam.MIPEmphasis, 3);
			// cplex.setParam(IloCplex.IntParam.pruning, 1);
			// formulation1.setParam(IloCplex.Param.MIP.Strategy.File,3); // not
			// needed when Emphasis.Memory==true
			// cplex.setParam(IloCplex.Param.Emphasis.MIP, 3);
			cplex.setParam(IloCplex.Param.WorkMem, 1024);

			cplex.solve();
			System.out.println("The obj= " + cplex.getObjValue());

			// System.out.println(cplex.solve());

		} catch (IloException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public void output() throws FileNotFoundException, UnknownObjectException, IloException {
		// output the solution
		PrintWriter out = new PrintWriter("outcome_check.txt");

		// X[i][j][k]
		// out.println("X[i][j][k]");
		out.println("Variable X");
		// for (int i = 0; i < numberOfCities * (T + 3); i++) {
		// for (int j = 0; j < numberOfCities * (T + 3); j++) {
		// for (int k = 0; k < numberOfTrucks; k++) {
		// if (connect[i][j]) {
		// if (Math.abs(cplex.getValue(X[i][j][k]) - 0) > e) {
		// out.println("vehicle " + k + ": " + getName(i) + "->" + getName(j) + " = "
		// + cplex.getValue(X[i][j][k]));
		// }
		// // out.print(cplex.getValue(X[i][j][k]) + " ");
		// }
		// }
		// }
		// }

		for (int k = 0; k < numberOfTrucks; k++) {
			for (int i = 0; i < numberOfCities * (T + 3); i++) {
				for (int j = 0; j < numberOfCities * (T + 3); j++) {
					if (connect[i][j]) {
						if (Math.abs(cplex.getValue(X[i][j][k]) - 0) > e) {
							out.println("vehicle " + k + ": " + getName(i) + "->" + getName(j) + " = "
									+ cplex.getValue(X[i][j][k]));
						}
					}
				}
			}
		}

		// y[i][j][p]
		// out.println();
		out.println("y[i][j][p]");
		// for (int i = 0; i < numberOfCities * (T + 1); i++) {
		// for (int j = 0; j < numberOfCities * (T + 1); j++) {
		//
		// for (int p = 0; p < numberOfDemandPair; p++) {
		// if (connect[i][j]) {
		// if (Math.abs(cplex.getValue(y[i][j][p]) - 0) > e) {
		// out.println(getName(i) + "->" + getName(j) + "for demand " +
		// demandPairs.get(p).s + "->"
		// + demandPairs.get(p).t + " :" + cplex.getValue(y[i][j][p]));
		//
		// // out.print(cplex.getValue(x[i][j][k][od]) + " ");
		// }
		// }
		// }
		//
		// }
		// }

		for (int p = 0; p < numberOfDemandPair; p++) {
			out.println("for demand " + demandPairs.get(p).s + "->" + demandPairs.get(p).t);
			for (int j = 0; j < numberOfCities * (T + 1); j++) {

				for (int i = 0; i < numberOfCities * (T + 1); i++) {
					if (connect[i][j]) {
						if (Math.abs(cplex.getValue(y[i][j][p]) - 0) > e) {
							out.println(getName(i) + "->" + getName(j) + " :" + cplex.getValue(y[i][j][p]));

							// out.print(cplex.getValue(x[i][j][k][od]) + " ");
						}
					}
				}

			}
		}

		out.close();
	}

	public static void main(String[] args) throws IOException, UnknownObjectException, IloException {
		MulticommodityFlowModel test = new MulticommodityFlowModel();
		// test.readData("out_small.txt");
		// test.readData("./data/temp.txt");
		// test.readData("./data/out_small3.txt");
		// test.readData("./data/out_small2.txt");
		// test.readData("./data/data2.txt");
		// test.readData("./data/data1_1.txt");
		// test.readData("./data/out2.txt");
		// test.readData("./data/out_small3_4.txt");
		// test.readData("./data/check10_50_3.txt");
		test.readData("./data/report4_4.txt");
		test.graphTransfer();
		test.ModelBuilding();
		test.output();
	}
}
