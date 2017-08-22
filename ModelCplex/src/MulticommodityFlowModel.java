import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;



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
	
	private ArrayList<ArrayList<Integer>> distance;
	
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
 * 
 */
	public void graphTransfer() {
		
	}
	
	public static void main(String[] args) throws IOException {
		MulticommodityFlowModel test=new MulticommodityFlowModel();
		test.readData("out2.txt");
	}
}
