import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashSet;

public class DataGeneration {
    //city parameter
    private int numberOfCity=4;
    private int numberOfDemandODs=3;
    private int demandLb=2000;
    private int demandUb=5000;
    private int coordinateLimitx=1200;
    private int coordinateLimity=1200;
    private int openTime=8;
    private int closeTime=20;
    private int arrivalTime=4;
    private double processingTime=1;
    private int M=1000000;
    
    //truck parameter
    private double fixedCost=2000;
    private double transportationCost=0.002;
    private int legLimit=2;
    private double distanceLimit=5000;
    private double speed=6000;
    private double drivingTimeLimit=10;
    private int numberOfTrucksPerCity=1;
    private int truckCapacity=5000;
    
    
    
    private PrintWriter out;
    
    public DataGeneration(String out) throws FileNotFoundException{
        this.out=new PrintWriter(out);
    }
    
    public void close(){
        out.close();
    }
    
    public void NodeParameter(){
        // declare all the citys
        out.print("Citys={");
        for(int i=0;i<numberOfCity;i++){
            if(i<numberOfCity-1){
                out.print("'"+toStringCity(i)+"',");
            }else{
                out.print("'"+toStringCity(i)+"'");
            }
        }
        out.println("}");
        
        
        //declare all the demandODs and demandQuantity
        out.println("demandQuantity = {");
        HashSet<Integer> set=new HashSet<Integer>();
        for(int i=0;i<numberOfDemandODs;i++){
            int s=(int) (Math.random()*numberOfCity);
            int t=(int) (Math.random()*numberOfCity);
            int check=s*1000+t;
            
            while(set.contains(check)||s==t){
                s=(int) (Math.random()*numberOfCity);
                t=(int) (Math.random()*numberOfCity);
                check=s*1000+t;
            }
            
            int quantity=(int) (Math.random()*(demandUb-demandLb)+demandLb);
            out.println("('"+toStringCity(s)+"','"+toStringCity(t)+"'): "+quantity);
            set.add(check);
            
        }
        out.println("}");
        
        //declare the coordinate of cities
        out.println("coordinate={");
        for(int i=0;i<numberOfCity;i++){
            double x=Math.random()*coordinateLimitx;
            double y=Math.random()*coordinateLimity;
            out.println("'"+toStringCity(i)+"':("+x+","+y+")");
        }
        out.println("}");
        
        //open time of all cities(HLi)
        out.print("openTime={");
        for(int i=0;i<numberOfCity;i++){
            if(i<numberOfCity-1){
                out.print("('"+toStringCity(i)+"'):"+openTime+",");
            }else{
                out.print("('"+toStringCity(i)+"'):"+openTime);
            }

        }
        out.println("}");
        
        //close time of all cities(HUi)
        out.print("closeTime={");
        for(int i=0;i<numberOfCity;i++){
            if(i<numberOfCity-1){
                out.print("('"+toStringCity(i)+"'):"+closeTime+",");
            }else{
                out.print("('"+toStringCity(i)+"'):"+closeTime);
            }

        }
        out.println("}");
        
        //arrival time lower bound of all cities(ai)
        out.print("arrivalTimeLB={");
        for(int i=0;i<numberOfCity;i++){
            if(i<numberOfCity-1){
                out.print("('"+toStringCity(i)+"'):"+arrivalTime+",");
            }else{
                out.print("('"+toStringCity(i)+"'):"+arrivalTime);
            }

        }
        out.println("}");
        
        //processing time of all cities(i)
        out.print("processingTime={");
        for(int i=0;i<numberOfCity;i++){
            if(i<numberOfCity-1){
                out.print("('"+toStringCity(i)+"'):"+processingTime+",");
            }else{
                out.print("('"+toStringCity(i)+"'):"+processingTime);
            }

        }
        out.println("}");
        
        //sufficiently large value M
        out.println("M="+M);
        
        out.println();
        
    }
    
    public void TruckParameter(){
        
        out.println("fixedCost="+fixedCost);
        out.println("transportationCost="+transportationCost);
        out.println("max number of legs per truck="+legLimit);
        out.println("max distance by a truck="+distanceLimit);
        out.println("averageSpeed="+speed);
        out.println("driving time per day="+drivingTimeLimit);
        
        //trucks
        out.println("trucks={");
        for(int i=0;i<numberOfCity*numberOfTrucksPerCity;i++){
            if(i%numberOfTrucksPerCity==numberOfTrucksPerCity-1){
                out.println("'t"+toStringCity(i/numberOfTrucksPerCity)+"-"+i%numberOfTrucksPerCity+"',");
            }else{
                out.print("'t"+toStringCity(i/numberOfTrucksPerCity)+"-"+i%numberOfTrucksPerCity+"',");
            }
        }
        out.println("}");
        
        //truckCapacity
        out.println("truckCapacity={");
        for(int i=0;i<numberOfCity*numberOfTrucksPerCity;i++){
            if(i%numberOfTrucksPerCity==numberOfTrucksPerCity-1){
                out.println("'t"+toStringCity(i/numberOfTrucksPerCity)+"-"+i%numberOfTrucksPerCity+"':"+truckCapacity+",");
            }else{
                out.print("'t"+toStringCity(i/numberOfTrucksPerCity)+"-"+i%numberOfTrucksPerCity+"':"+truckCapacity+",");
            }
        }
        out.println("}");
        
        //truck startNode
        out.println("truckStartNode={");
        for(int i=0;i<numberOfCity*numberOfTrucksPerCity;i++){
            if(i%numberOfTrucksPerCity==numberOfTrucksPerCity-1){
                out.println("'t"+toStringCity(i/numberOfTrucksPerCity)+"-"+i%numberOfTrucksPerCity+"':"+toStringCity(i/numberOfTrucksPerCity)+",");
            }else{
                out.print("'t"+toStringCity(i/numberOfTrucksPerCity)+"-"+i%numberOfTrucksPerCity+"':"+toStringCity(i/numberOfTrucksPerCity)+",");
            }
        }
        out.println("}");
        
        
        
    }
    
    public String toStringCity(int city){
        String string="";
        if(city<10){
            string="00"+city;
        }else{
            if(city>=10&&city<=99){
                string="0"+city;
            }else{
                string=String.valueOf(city);
            }
        }
        
        return string;
    }

    public static void main(String[] args) throws FileNotFoundException {
        DataGeneration data=new DataGeneration("out_small.txt");
        data.NodeParameter();
        data.TruckParameter();
        data.close();
    }
}
