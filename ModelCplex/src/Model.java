import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Paths;
import java.util.*;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.*;

public class Model {
    private HashMap<String, Integer> cityIndex;
    private HashMap<String, Integer> truckIndex;
    private int numberOfCities;
    private int numberOfTrucks;
    // private HashSet<demandPair> demandPairs;
    private ArrayList<demandPair> demandPairs;
    private int numberOfDemandPair;
    private double x[];
    private double y[];
    private double[][] distance;
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
    // private ArrayList<Integer>[] nodeTrucks;
    private int[] truckStartNode;

    private IloCplex cplex;

    private class demandPair {
        private int s;
        private int t;
        private int demandQuantity;

    }

    public Model() throws IOException {
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

        // System.out.println(temp.toString());
        // int a=temp.lastIndexOf(": ");
        // System.out.println(a);
        // int demand=Integer.parseInt(temp.substring(a+2,temp.length()));
        // System.out.println(demand);

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
        x = new double[numberOfCities];
        y = new double[numberOfCities];
        distance = new double[numberOfCities][numberOfCities];
        temp = in.nextLine();

        while (!temp.equals("}")) {
            lIndex = temp.indexOf("'", 0);
            rIndex = temp.indexOf("'", lIndex + 1);
            int city = Integer.parseInt(temp.substring(lIndex + 1, rIndex));

            lIndex = temp.indexOf("(", rIndex + 1);
            rIndex = temp.indexOf(",", lIndex + 1);
            x[city] = Double.valueOf(temp.substring(lIndex + 1, rIndex));
            y[city] = Double.valueOf(temp.substring(rIndex + 1, temp.length() - 1));
            // System.out.println(city);
            // System.out.println("x="+x[city]);
            // System.out.println("y="+y[city]);

            temp = in.nextLine();
        }

        // Calculate distance
        for (int i = 0; i < numberOfCities; i++) {
            for (int j = i; j < numberOfCities; j++) {
                if (i == j) {
                    distance[i][j] = 0;
                } else {
                    distance[i][j] = Math.sqrt(Math.pow((x[i] - x[j]), 2) + Math.pow((y[i] - y[j]), 2));
                    distance[j][i] = distance[i][j];
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

        // ---read "closeTime"---//
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

        // ---read "closeTime"---//
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

        // System.out.println(Arrays.toString(truckStartNode));

        // private HashMap<String, Integer> cityIndex;
        // private HashMap<String, Integer> truckIndex;
        // private int numberOfCities;
        // private int numberOfTrucks;
        // private HashSet<demandPair> demandPairs;
        // private double x[];
        // private double y[];
        // private double[][] distance;
        // private double[] openTime;
        // private double[] closeTime;
        // private double[] arrivalTime;
        // private double[] processingTime;
        // private int M;
        // private double fixedCost;
        // private double transportationCost;
        // private int legLimit;
        // private double distanceLimit;
        // private double averageSpeed;
        // private double drivingTimePerDay;
        // private int[] truckCapacity;
        // private ArrayList<Integer>[] nodeTrucks;

        // ---Output---//
        // System.out.println("number of cities= "+numberOfCities);
        // System.out.println("number of trucks= "+numberOfTrucks);

    }

    public void ModelBuilding() throws FileNotFoundException {
        try {
            cplex = new IloCplex();
            // ---decision variables---//
            IloIntVar[][][] X = new IloIntVar[numberOfCities][numberOfCities][numberOfTrucks];
            IloIntVar[][][][] x = new IloIntVar[numberOfCities][numberOfCities][numberOfTrucks][numberOfDemandPair];
            IloNumVar[][][][] y = new IloNumVar[numberOfCities][numberOfCities][numberOfTrucks][numberOfDemandPair];
            IloNumVar[][] ta = new IloNumVar[numberOfTrucks][numberOfCities];
            IloNumVar[][] td = new IloNumVar[numberOfTrucks][numberOfCities];

            for (int i = 0; i < numberOfCities; i++) {
                for (int j = 0; j < numberOfCities; j++) {
                    for (int k = 0; k < numberOfTrucks; k++) {

                        if (i != j) {
                            X[i][j][k] = cplex.intVar(0, 1, "X" + i + "," + j + "," + k);
                        } else {
                            X[i][j][k] = cplex.intVar(0, 0, "X" + i + "," + j + "," + k);
                        }

                        // for (int s = 0; s < numberOfCities; s++) {
                        // for (int t = 0; t < numberOfCities; t++) {
                        // if (s == t || i == j) {
                        // x[i][j][k][s][t] = cplex.intVar(0, 0,
                        // "x" + i + "," + j + "," + k + "," + s + "," + t);
                        // y[i][j][k][s][t] = cplex.numVar(0, 0,
                        // "y" + i + "," + j + "," + k + "," + s + "," + t);
                        // } else {
                        // x[i][j][k][s][t] = cplex.intVar(0, 1,
                        // "x" + i + "," + j + "," + k + "," + s + "," + t);
                        // y[i][j][k][s][t] = cplex.numVar(0, Double.MAX_VALUE,
                        // "y" + i + "," + j + "," + k + "," + s + "," + t);
                        // }
                        // }
                        // }
                        for (int od = 0; od < numberOfDemandPair; od++) {
                            int s = demandPairs.get(od).s;
                            int t = demandPairs.get(od).t;
                            if (s == t || i == j) {
                                x[i][j][k][od] = cplex.intVar(0, 0, "x" + i + "," + j + "," + k + "," + s + "," + t);
                                y[i][j][k][od] = cplex.numVar(0, 0, "y" + i + "," + j + "," + k + "," + s + "," + t);
                            } else {
                                x[i][j][k][od] = cplex.intVar(0, 1, "x" + i + "," + j + "," + k + "," + s + "," + t);
                                y[i][j][k][od] = cplex.numVar(0, Double.MAX_VALUE,
                                        "y" + i + "," + j + "," + k + "," + s + "," + t);
                            }
                        }

                    }
                }
            }

            for (int k = 0; k < numberOfTrucks; k++) {
                for (int i = 0; i < numberOfCities; i++) {
                    ta[k][i] = cplex.numVar(0, Double.MAX_VALUE, "ta" + k + "," + i);
                    td[k][i] = cplex.numVar(0, Double.MAX_VALUE, "td" + k + "," + i);
                }
            }

            // ---Objective---//
            IloLinearNumExpr obj = cplex.linearNumExpr();

            double parameter = fixedCost / (averageSpeed * drivingTimePerDay);
            for (int i = 0; i < numberOfCities; i++) {
                for (int j = 0; j < numberOfCities; j++) {
                    for (int k = 0; k < numberOfTrucks; k++) {

                        obj.addTerm(parameter * distance[i][j], X[i][j][k]);

                    }
                }
            }

            for (int k = 0; k < numberOfTrucks; k++) {
                for (int od = 0; od < numberOfDemandPair; od++) {
                    for (int i = 0; i < numberOfCities; i++) {
                        for (int j = 0; j < numberOfCities; j++) {

                            obj.addTerm(transportationCost * distance[i][j], y[i][j][k][od]);

                        }
                    }
                }
            }
            cplex.addMinimize(obj);

            // ---constraint 1 & 3 & 4---//
            for (int k = 0; k < numberOfTrucks; k++) {

                IloLinearNumExpr constraint3 = cplex.linearNumExpr();
                IloLinearNumExpr constraint4 = cplex.linearNumExpr();

                for (int i = 0; i < numberOfCities; i++) {

                    IloLinearNumExpr constraint1 = cplex.linearNumExpr();

                    for (int j = 0; j < numberOfCities; j++) {
                        constraint1.addTerm(1, X[j][i][k]);
                        constraint1.addTerm(-1, X[i][j][k]);

                        constraint3.addTerm(1, X[i][j][k]);
                        constraint4.addTerm(distance[i][j], X[i][j][k]);
                    }
                    cplex.addEq(0, constraint1);
                }

                cplex.addGe(legLimit, constraint3);
                // System.out.println("Constrint3 is "+constraint3.toString());
                cplex.addGe(distanceLimit, constraint4);
            }

            // ---constraint 2---//
            for (int k = 0; k < numberOfTrucks; k++) {
                int startNode = truckStartNode[k];
                IloLinearNumExpr constraint2 = cplex.linearNumExpr();

                for (int j = 0; j < numberOfCities; j++) {
                    constraint2.addTerm(1, X[startNode][j][k]);
                }
                cplex.addGe(1, constraint2);
                // System.out.println(constraint2.toString());
            }

            // ---constraint 7-1,7-2,8-1,8-2 ---//
            // for(int k=0;k<numberOfTrucks;k++){
            // for(int i=0;i<numberOfCities;i++){
            // IloLinearNumExpr constraint7_1 = cplex.linearNumExpr();
            // IloLinearNumExpr constraint7_2 = cplex.linearNumExpr();
            // IloLinearNumExpr constraint8_1 = cplex.linearNumExpr();
            // IloLinearNumExpr constraint8_2 = cplex.linearNumExpr();
            // for(int j=0;j<numberOfCities;j++){
            // constraint7_1.addTerm((-1)*M, X[j][i][k]);
            // constraint7_2.addTerm((-1)*M, X[j][i][k]);
            // constraint8_1.addTerm(M, X[j][i][k]);
            // constraint8_2.addTerm(M, X[j][i][k]);
            // }
            //
            //
            // constraint7_1.addTerm(-1,td[k][i]%24);
            // constraint7_2.addTerm(-1,ta[k][i]%24);
            // constraint8_1.addTerm(-1,td[k][i]);
            // constraint8_2.addTerm(-1,ta[k][i]);
            //
            // cplex.addLe(-closeTime[], e)
            // }
            // }

            // ---constraint 5 ---//
            for (int k = 0; k < numberOfTrucks; k++) {
                for (int i = 0; i < numberOfCities; i++) {
                    IloLinearNumExpr constraint5 = cplex.linearNumExpr();
                    for (int j = 0; j < numberOfCities; j++) {
                        constraint5.addTerm(M, X[j][i][k]);
                    }
                    constraint5.addTerm(1, td[k][i]);
                    // System.out.println(constraint5.toString());
                    cplex.addGe(arrivalTime[i] + M, constraint5);
                }
            }

            // ---constraint 6 (lack)---//

            // ---constraint 9 & 10 ---//
            for (int k = 0; k < numberOfTrucks; k++) {
                for (int i = 0; i < numberOfCities; i++) {
                    for (int j = 0; j < numberOfCities; j++) {
                        if (i != j) {
                            IloLinearNumExpr constraint9 = cplex.linearNumExpr();
                            IloLinearNumExpr constraint10 = cplex.linearNumExpr();

                            constraint9.addTerm(1, ta[k][j]);
                            constraint10.addTerm(1, ta[k][j]);

                            constraint9.addTerm(-1, td[k][i]);
                            constraint10.addTerm(-1, td[k][i]);

                            constraint9.addTerm(M, X[i][j][k]);
                            constraint10.addTerm((-1) * M, X[i][j][k]);

                            cplex.addGe((distance[i][j] / averageSpeed) + M, constraint9);
                            cplex.addLe((distance[i][j] / averageSpeed) - M, constraint10);

                            // System.out.println(constraint9.toString());
                            // System.out.println(constraint10.toString());
                        }
                    }
                }
            }

            // ---constraint 11 ---//
            for (int k = 0; k < numberOfTrucks; k++) {
                for (int i = 0; i < numberOfCities; i++) {
                    if (i != truckStartNode[k]) {
                        IloLinearNumExpr constraint11 = cplex.linearNumExpr();
                        constraint11.addTerm(1, td[k][i]);
                        constraint11.addTerm(-1, ta[k][i]);
                        cplex.addEq(processingTime[i], constraint11);
                        // System.out.println(constraint11.toString());
                    }
                }
            }

            // ---constraint 12 ---//
            for (int i = 0; i < numberOfCities; i++) {
                for (int od = 0; od < numberOfDemandPair; od++) {
                    for (int k1 = 0; k1 < numberOfTrucks; k1++) {
                        for (int k2 = 0; k2 < numberOfTrucks; k2++) {
                            if (k1 != k2) {
                                for (int i1 = 0; i1 < numberOfCities; i1++) {
                                    for (int i2 = 0; i2 < numberOfCities; i2++) {
                                        if (i1 != i && i2 != i) {
                                            IloLinearNumExpr constraint12 = cplex.linearNumExpr();
                                            constraint12.addTerm(1, ta[k1][i]);
                                            constraint12.addTerm(-1, td[k2][i]);
                                            constraint12.addTerm(M, x[i1][i][k1][od]);
                                            constraint12.addTerm(M, x[i][i2][k2][od]);

                                            cplex.addGe(2 * M - processingTime[i], constraint12);
                                            // System.out.println(constraint12);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // ---constraint 13 & 14 ---//
            for (int od = 0; od < numberOfDemandPair; od++) {
                int s = demandPairs.get(od).s;
                int t = demandPairs.get(od).t;

                IloLinearNumExpr constraint13 = cplex.linearNumExpr();
                IloLinearNumExpr constraint14 = cplex.linearNumExpr();

                for (int k = 0; k < numberOfTrucks; k++) {
                    for (int j = 0; j < numberOfCities; j++) {
                        constraint13.addTerm(1, y[s][j][k][od]);
                        constraint14.addTerm(1, y[j][t][k][od]);
                    }
                }

                cplex.addEq(demandPairs.get(od).demandQuantity, constraint13);
                cplex.addEq(demandPairs.get(od).demandQuantity, constraint14);
                // System.out.println(constraint14.toString());
            }

            // ---constraint 15 ---//
            for (int od = 0; od < numberOfDemandPair; od++) {

                int s = demandPairs.get(od).s;
                int t = demandPairs.get(od).t;

                for (int i = 0; i < numberOfCities; i++) {
                    if (i != s && i != t) {
                        IloLinearNumExpr constraint15 = cplex.linearNumExpr();
                        for (int k = 0; k < numberOfTrucks; k++) {
                            for (int j = 0; j < numberOfCities; j++) {
                                constraint15.addTerm(1, y[i][j][k][od]);
                                constraint15.addTerm(-1, y[j][i][k][od]);
                            }
                        }
                        cplex.addEq(0, constraint15);
                        // System.out.println(constraint15.toString());
                    }
                }
            }

            // ---constraint 16 & 17 ---//
            for (int k = 0; k < numberOfTrucks; k++) {
                for (int od = 0; od < numberOfDemandPair; od++) {
                    for (int i = 0; i < numberOfCities; i++) {
                        for (int j = 0; j < numberOfCities; j++) {
                            if (i != j) {
                                IloLinearNumExpr constraint16 = cplex.linearNumExpr();
                                IloLinearNumExpr constraint17 = cplex.linearNumExpr();

                                constraint16.addTerm(1, y[i][j][k][od]);
                                constraint16.addTerm(-M, x[i][j][k][od]);
                                cplex.addGe(0, constraint16);
                                // System.out.println(constraint16.toString());

                                constraint17.addTerm(1, x[i][j][k][od]);
                                constraint17.addTerm(-1, X[i][j][k]);
                                cplex.addGe(0, constraint17);
                                // System.out.println(constraint17.toString());
                            }
                        }
                    }
                }
            }

            // ---constraint 18 ---//
            for (int k = 0; k < numberOfTrucks; k++) {
                for (int od = 0; od < numberOfDemandPair; od++) {
                    IloLinearNumExpr constraint18 = cplex.linearNumExpr();
                    int s = demandPairs.get(od).s;
                    for (int j = 0; j < numberOfCities; j++) {
                        constraint18.addTerm(1, x[s][j][k][od]);
                    }

                    cplex.addGe(1, constraint18);
                    // System.out.println(constraint18.toString());
                }
            }

            // ---constraint 19 ---//
            for (int k = 0; k < numberOfTrucks; k++) {
                for (int i = 0; i < numberOfCities; i++) {
                    for (int j = 0; j < numberOfCities; j++) {
                        IloLinearNumExpr constraint19 = cplex.linearNumExpr();

                        for (int od = 0; od < numberOfDemandPair; od++) {
                            constraint19.addTerm(1, y[i][j][k][od]);
                        }

                        cplex.addGe(truckCapacity[k], constraint19);
                        // System.out.println(constraint19.toString());

                    }
                }
            }

            cplex.solve();
            // System.out.println(cplex.solve());
            // cplex.exportModel("test.lp");
            // System.out.println(cplex.toString());

            // output the solution
            PrintWriter out = new PrintWriter("outcome.txt");

            for (int i = 0; i < numberOfCities; i++) {
                for (int j = 0; j < numberOfCities; j++) {
                    for (int k = 0; k < numberOfTrucks; k++) {
                        if (cplex.getValue(X[i][j][k]) != 0) {
                            out.println(X[i][j][k].getName() + "=" + cplex.getValue(X[i][j][k]));
//                            System.out.println(X[i][j][k].getName() + "=" + cplex.getValue(X[i][j][k]));
                        }
                    }
                }
            }

            // for (int i = 0; i < numberOfCities; i++) {
            // for (int j = 0; j < numberOfCities; j++) {
            // if(i!=j){
            // for (int k = 0; k < numberOfTrucks; k++) {
            // for (int od = 0; od < numberOfDemandPair; od++) {
            //// System.out.println("i=" + i + ", j=" + j + ", k=" + k + ", od="
            // + od);
            // if (cplex.getValue(x[i][j][k][od]) != 0){
            // out.println(x[i][j][k][od].getName() + "=" +
            // cplex.getValue(x[i][j][k][od]));
            //// System.out.println(x[i][j][k][od].getName() + "=" +
            // cplex.getValue(x[i][j][k][od]));
            // }
            // }
            // }
            // }
            //
            // }
            // }
            for (int i = 0; i < numberOfCities; i++) {
                for (int j = 0; j < numberOfCities; j++) {

                    for (int k = 0; k < numberOfTrucks; k++) {
                        for (int od = 0; od < numberOfDemandPair; od++) {
                            System.out.println("i=" + i + ", j=" + j + ", k=" + k + ", od=" + od);
                            if (cplex.getValue(x[i][j][k][od]) != 0) {
                                out.println(x[i][j][k][od].getName() + "=" + cplex.getValue(x[i][j][k][od]));
                                // System.out.println(x[i][j][k][od].getName() +
                                // "=" + cplex.getValue(x[i][j][k][od]));
                            }
                        }
                    }

                }
            }

            for (int i = 0; i < numberOfCities; i++) {
                for (int j = 0; j < numberOfCities; j++) {
                    for (int k = 0; k < numberOfTrucks; k++) {
                        for (int od = 0; od < numberOfDemandPair; od++) {
                            if (cplex.getValue(y[i][j][k][od]) != 0)
                                out.println(y[i][j][k][od].getName() + "=" + cplex.getValue(y[i][j][k][od]));
                        }
                    }
                }
            }

            for (int k = 0; k < numberOfTrucks; k++) {
                for (int i = 0; i < numberOfCities; i++) {
                    if (cplex.getValue(ta[k][i]) != 0)
                        out.println(ta[k][i].getName() + "=" + cplex.getValue(ta[k][i]));
                }
            }

            for (int k = 0; k < numberOfTrucks; k++) {
                for (int i = 0; i < numberOfCities; i++) {
                    if (cplex.getValue(td[k][i]) != 0)
                        out.println(td[k][i].getName() + "=" + cplex.getValue(td[k][i]));
                }
            }

            out.close();

        } catch (IloException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    public static void main(String[] args) throws IOException {
        Model model = new Model();
        model.readData("out.txt");
        model.ModelBuilding();
    }
}
