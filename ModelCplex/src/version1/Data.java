package version1;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.*;




public class Data {
	
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
    
    private final int T = 36;
    private int[][] b;
    static double e=Math.pow(10, -2);
    private boolean[][] connect;

    private class Edge {
        int start,end;
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

    private ArrayList<Edge> edgeSet;
    private ArrayList<ArrayList<Integer>> distance;
    private ArrayList<ArrayList<Integer>> distanceReverse;
    private int numOfEdge1,numOfEdge2,numOfEdge3,numOfEdge4;
    
    
    private double[] c,f,bb;
    private double[][] A,B;
    
    public Data() {
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

        edgeSet=new ArrayList<Edge>();
        distance = new ArrayList<ArrayList<Integer>>();
        distanceReverse = new ArrayList<ArrayList<Integer>>();

        // only record Vst and O
        for (int i = 0; i < numberOfCities * (T + 2); i++) {
            ArrayList<Integer> templist = new ArrayList<Integer>();
            distance.add(templist);

        }

        for (int i = 0; i < numberOfCities * (T + 3); i++) {
            ArrayList<Integer> templist2 = new ArrayList<Integer>();
            distanceReverse.add(templist2);
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
                        edge.start=nodeIndex1;
                        edge.end = nodeIndex2;
                        edge.length = length[i][j];
                        edge.setIndex = 1;
                        edge.u = i;
                        edge.v = j;
                        edge.t1 = t - timeLength;
                        edge.t2 = t;
                        edgeSet.add(edge);
                        distance.get(nodeIndex1).add(edgeSet.size()-1);
                        distanceReverse.get(nodeIndex2).add(edgeSet.size()-1);
                        connect[nodeIndex1][nodeIndex2] = true;

                        t++;
                        nodeIndex1++;
                        nodeIndex2++;
                    }
                }
            }
        }
        numOfEdge1=edgeSet.size();

        // add AT
        for (int node = 0; node < numberOfCities; node++) {
            for (int t = 0; t < T; t++) {
                Edge edge = new Edge();
                int nodeIndex = node * (T + 1) + t;
                edge.start=nodeIndex;
                edge.end= nodeIndex + 1;
                edge.length = 0;
                edge.setIndex = 2;
                edge.u = node;
                edge.v = node;
                edge.t1 = t;
                edge.t2 = t + 1;
                edgeSet.add(edge);
                
                distance.get(nodeIndex).add(edgeSet.size()-1);
                distanceReverse.get(nodeIndex + 1).add(edgeSet.size()-1);
                connect[nodeIndex][nodeIndex + 1] = true;
            }
        }
        
        numOfEdge2=edgeSet.size()-numOfEdge1;


        // add AO:(Ok,n0)
        for (int o = 0; o < numberOfCities; o++) {
            int oIndex = numberOfCities * (T + 1) + o;

            for (int node = 0; node < numberOfCities; node++) {
                int nodeIndex = node * (T + 1);
                Edge edge = new Edge();
                edge.start=oIndex;
                edge.end = nodeIndex;
                edge.length = length[o][node];
                edge.setIndex = 3;
                edge.u = o;
                edge.v = node;
                edge.t1 = -1;
                edge.t2 = 0;
                
                edgeSet.add(edge);
                distance.get(oIndex).add(edgeSet.size()-1);
                distanceReverse.get(edge.end).add(edgeSet.size()-1);
                connect[oIndex][nodeIndex] = true;
            }
        }
        
        numOfEdge3=edgeSet.size()-numOfEdge1-numOfEdge2;

        // add AD:(nT,Dk)
        for (int node = 0; node < numberOfCities; node++) {
            int nodeIndex = node * (T + 1) + T;
            for (int d = 0; d < numberOfCities; d++) {
                int dIndex = numberOfCities * (T + 2) + d;
                Edge edge = new Edge();
                edge.start=nodeIndex;
                edge.end = dIndex;
                edge.length = length[node][d];
                edge.setIndex = 4;
                edge.u = node;
                edge.v = d;
                edge.t1 = T;
                edge.t2 = -1;
                edgeSet.add(edge);
                distance.get(nodeIndex).add(edgeSet.size()-1);
                distanceReverse.get(dIndex).add(edgeSet.size()-1);
                connect[nodeIndex][dIndex] = true;

            }
        }
        
        numOfEdge4=edgeSet.size()-numOfEdge1-numOfEdge2-numOfEdge3;

        // set b
        b = new int[numberOfCities * (T + 1)][numberOfDemandPair];
        for (int p = 0; p < numberOfDemandPair; p++) {
            demandPair pair = demandPairs.get(p);
            b[pair.s * (T + 1)][p] = pair.demandQuantity;
            b[pair.t * (T + 1) + T][p] = -pair.demandQuantity;
        }


    }
	
    public void matrixGenerator() {
    	int numOfEdge12=numOfEdge1+numOfEdge2;
    	int numOfy=(numOfEdge12)*numberOfDemandPair;
    	int numOfx=edgeSet.size()*numberOfTrucks;
    	
    	
    	c=new double[numOfy];
    	int index=0;
    	for(int p=0;p<numberOfDemandPair;p++) {
    		for(int e=0;e<numOfEdge12;e++) {
    			c[index]=transportationCost*edgeSet.get(e).length;
    			index++;
    		}
    	}
    	
//    	System.out.println(Arrays.toString(c));
    	
    	f=new double[numOfx];
    	double constant=fixedCost/(averageSpeed*drivingTimePerDay);
        index=0;
    	for(int k=0;k<numberOfTrucks;k++) {
    		for(int e=0;e<edgeSet.size();e++) {
    			f[index]=constant*edgeSet.get(e).length;
    			index++;
    		}
    	}
    	
//    	System.out.println(Arrays.toString(f));
//    	System.out.println(numOfEdge1);
//    	System.out.println(numOfy);
    	
    	
    	A=new double[numberOfCities*(T+1)*numberOfDemandPair*2+numOfEdge1][numOfy];
    	B=new double[numberOfCities*(T+1)*numberOfDemandPair*2+numOfEdge1][numOfx];
    	bb=new double[numberOfCities*(T+1)*numberOfDemandPair*2+numOfEdge1];
    	
    	
    	int row=0;
    	for(int p=0;p<numberOfDemandPair;p++) {
    		for(int i=0;i<numberOfCities*(T+1);i++) {
    			
    			
//    			for(int e=0;e<numOfEdge12;e++) {
//    				if(edgeSet.get(e).start==i) {
//    					A[row][numOfEdge12*p+e]=1;
//    				}else if(edgeSet.get(e).end==i) {
//    					A[row][numOfEdge12*p+e]=-1;
//    				}
//    			}
    			
    			for(int e=0;e<distance.get(i).size();e++) {
    				if(edgeSet.get(distance.get(i).get(e)).setIndex<3) {
    					A[row][numOfEdge12*p+distance.get(i).get(e)]=1;
    				}
    			}
    			
    			for(int e=0;e<distanceReverse.get(i).size();e++) {
    				if(edgeSet.get(distanceReverse.get(i).get(e)).setIndex<3) {
    					A[row][numOfEdge12*p+distanceReverse.get(i).get(e)]=-1;
    				}
    			}
    			
    			
    			bb[row]=b[i][p];
    			row++;
    			
    			
    		}
    	}
    	
//    	System.out.println("node 0 point to: ");
//    	for(int e=0;e<edgeSet.size();e++) {
//    		if(edgeSet.get(e).start==0) {
//    			System.out.print(e+" ");
//    		}
//    	}
//    	
//    	System.out.println();
//    	System.out.println("node 0 point from: ");
//    	for(int e=0;e<edgeSet.size();e++) {
//    		if(edgeSet.get(e).end==0) {
//    			System.out.print(edgeSet.get(e).start+" "+edgeSet.get(e).setIndex+" ");
//    		}
//    	}
//    	
//    	
//    	System.out.println();
//    	int temp=0;
//    	for(int p=0;p<numberOfDemandPair;p++) {
//    		for(int e=0;e<numOfEdge12;e++) {
//    			System.out.print(A[0][temp]+" ");
//    			temp++;
//    		}
//    		System.out.println();
//    	}
//    	System.out.println();
//    	temp=0;
//    	for(int p=0;p<numberOfDemandPair;p++) {
//    		for(int e=0;e<numOfEdge12;e++) {
//    			System.out.print(A[222][temp]+" ");
//    			temp++;
//    		}
//    		System.out.println();
//    	}
    	
    	
    	
    	
    	for(int p=0;p<numberOfDemandPair;p++) {
    		for(int i=0;i<numberOfCities*(T+1);i++) {
    			
    			
//    			for(int e=0;e<numOfEdge12;e++) {
//    				if(edgeSet.get(e).start==i) {
//    					A[row][numOfEdge12*p+e]=-1;
//    				}else if(edgeSet.get(e).end==i) {
//    					A[row][numOfEdge12*p+e]=1;
//    				}
//    			}
    			
    			for(int e=0;e<distance.get(i).size();e++) {
    				if(edgeSet.get(distance.get(i).get(e)).setIndex<3) {
    					A[row][numOfEdge12*p+distance.get(i).get(e)]=-1;
    				}
    			}
    			
    			for(int e=0;e<distanceReverse.get(i).size();e++) {
    				if(edgeSet.get(distanceReverse.get(i).get(e)).setIndex<3) {
    					A[row][numOfEdge12*p+distanceReverse.get(i).get(e)]=1;
    				}
    			}
    			
    			bb[row]=-b[i][p];
    			row++;
    			
//    			if(p==0&&i==0) {
//    				System.out.println(Arrays.toString(A[row-1]));
//    			}
    			
    		}
    		
    	}
    	
    	
    	
    	for(int e=0;e<numOfEdge1;e++) {
    		
    		for(int p=0;p<numberOfDemandPair;p++) {
    			A[row][p*numOfEdge12+e]=-1;
    		}
    		
    		for(int k=0;k<numberOfTrucks;k++) {
    			B[row][k*edgeSet.size()+e]=truckCapacity[k];
    		}
    		
//    		if(e==0) {
//            	int temp=0;
//            	for(int k=0;k<numberOfTrucks;k++) {
//            		for(int ee=0;ee<edgeSet.size();ee++) {
//            			System.out.print(B[row][temp]+" ");
//            			temp++;
//            		}
//            		System.out.println();
//            	}
//    		}

    		row++;
    	}
    	

    	
    	
    	
    	
    	
    	
    }
	public static void main(String[] args) throws IOException {
		Data data=new Data();
		data.readData("out2.txt");
		data.graphTransfer();
		data.matrixGenerator();
	}
}