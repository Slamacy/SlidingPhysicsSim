package CA2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;
import java.util.Scanner;

/**
 * @authors
 * Patricia Bere - D00193593
 * Oisin Murphy - D00191700
 */

public class Main {
    private static double initialTime;
    private static double currentTime;
    private static double finalTime;
    private static double h = 0.1;
    private static double mass;
    private static double area;
    private static double volume;
    private static Vector<Double> forceGravity = new Vector<>(),
                                  position = new Vector<>(),
                                  velocity = new Vector<>(),
                                  forceGn = new Vector<>(),
                                  forceGp = new Vector<>(),
                                  forceNormal = new Vector<>(),
                                  forceFriction = new Vector<>(),
                                  forceNet = new Vector<>(),
                                  acceleration = new Vector<>(),
                                  kHat = new Vector<>(),
                                  nHat = new Vector<>(),
                                  forceGPHat = new Vector<>(),
                                  nBar = new Vector<>();
    private static double muStatic;
    private static double muKinetic;
    private static double n;
    private static double fGdotN;
    private static double forceNormalLength;
    private static double forceGPLength;
    private static boolean isFirst = true;
    private static Scanner in = new Scanner(System.in);
    private static File output = new File("output.txt");
    private static final double GRAVITY = 9.81;
    
    
    public static void main(String[] args) {
        kHat.add(0, 0.0);
        kHat.add(1, 0.0);
        kHat.add(2, 1.0);
        
        initialize();
        
        for(double i = initialTime; i < finalTime; i += h) {
            currentTime = i;
            
            gravForces();
            normalForce();
            frictionForce();
            netForce();
            getAcceleration();
            writeStepToFile();
            updatePosition();
            updateVelocity();
        }
    }
    
    public static void gravForces() {
        //forceGravity part Can be copied from last CA
        double fg1 = -(mass * GRAVITY * kHat.get(0));
        double fg2 = -(mass * GRAVITY * kHat.get(1));
        double fg3 = -(mass * GRAVITY * kHat.get(2));
        forceGravity.add(0, fg1);
        forceGravity.add(1, fg2);
        forceGravity.add(2, fg3);
        
        //get forceGn
        n = Math.sqrt((nBar.get(0) * nBar.get(0)) + (nBar.get(1) * nBar.get(1)) + (nBar.get(2) * nBar.get(2)));
        double nH1 = (1 / n) * nBar.get(0);
        double nH2 = (1 / n) * nBar.get(1);
        double nH3 = (1 / n) * nBar.get(2);
        
        nHat.add(0, nH1);
        nHat.add(1, nH2);
        nHat.add(2, nH3);
        
        fGdotN = (forceGravity.get(0) * nHat.get(0)) + (forceGravity.get(1) * nHat.get(1)) + (forceGravity.get(2) * nHat.get(2));
        if(fGdotN < 0) {
            fGdotN *= -1;
            //always has to be positive
        }
        
        double fGn1 = (fGdotN * nHat.get(0));
        double fGn2 = (fGdotN * nHat.get(1));
        double fGn3 = (fGdotN * nHat.get(2));
        
        forceGn.add(0, fGn1);
        forceGn.add(1, fGn2);
        forceGn.add(2, fGn3);
        
        //get forceGp
        double fGp1 = forceGravity.get(0) - forceGn.get(0);
        double fGp2 = forceGravity.get(1) - forceGn.get(1);
        double fGp3 = forceGravity.get(2) - forceGn.get(2);
        
        forceGp.add(0, fGp1);
        forceGp.add(1, fGp2);
        forceGp.add(2, fGp3);
        
        forceGPLength = Math.sqrt((forceGp.get(0) * forceGp.get(0)) + (forceGp.get(1) * forceGp.get(1)) + (forceGp.get(2) * forceGp.get(2)));
        
        double fGPhat1 = (1 / forceGPLength) * forceGp.get(0);
        double fGPhat2 = (1 / forceGPLength) * forceGp.get(1);
        double fGPhat3 = (1 / forceGPLength) * forceGp.get(2);
        
        forceGPHat.add(0, fGPhat1);
        forceGPHat.add(1, fGPhat2);
        forceGPHat.add(2, fGPhat3);
    }
    
    public static void normalForce() {
        double fN1 = forceGn.get(0) * -1;
        double fN2 = forceGn.get(1) * -1;
        double fN3 = forceGn.get(2) * -1;
        
        forceNormal.add(0, fN1);
        forceNormal.add(1, fN2);
        forceNormal.add(2, fN3);
    }
    
    public static void frictionForce() {
        forceNormalLength = Math.sqrt((forceNormal.get(0) * forceNormal.get(0)) + (forceNormal.get(1) * forceNormal.get(1)) + (forceNormal.get(2) * forceNormal.get(2)));
        
        double fF1;
        double fF2;
        double fF3;
        
        if(muStatic * forceNormalLength > fGdotN) {
            fF1 = (muKinetic * forceNormalLength) * forceGPHat.get(0);
            fF2 = (muKinetic * forceNormalLength) * forceGPHat.get(1);
            fF3 = (muKinetic * forceNormalLength) * forceGPHat.get(2);
            
            forceFriction.add(0, fF1);
            forceFriction.add(1, fF2);
            forceFriction.add(2, fF3);
        } else {
            fF1 = (muStatic * forceNormalLength) * forceGPHat.get(0);
            fF2 = (muStatic * forceNormalLength) * forceGPHat.get(1);
            fF3 = (muStatic * forceNormalLength) * forceGPHat.get(2);
            
            forceFriction.add(0, fF1);
            forceFriction.add(1, fF2);
            forceFriction.add(2, fF3);
        }
    }
    
    public static void netForce() {
        //Fnet = Fgn + Fgp + Fn + Ff
        double Fnet0 = forceGp.get(0) + forceGn.get(0) + forceNormal.get(0) + forceFriction.get(0);
        double Fnet1 = forceGp.get(1) + forceGn.get(1) + forceNormal.get(1) + forceFriction.get(1);
        double Fnet2 = forceGp.get(2) + forceGn.get(2) + forceNormal.get(2) + forceFriction.get(2);
        
        forceNet.add(0, Fnet0);
        forceNet.add(1, Fnet1);
        forceNet.add(2, Fnet2);
    }
    
    public static void getAcceleration() {
        double aHat0 = forceNet.get(0) / mass;
        double aHat1 = forceNet.get(1) / mass;
        double aHat2 = forceNet.get(2) / mass;
        
        acceleration.add(0, aHat0);
        acceleration.add(1, aHat1);
        acceleration.add(2, aHat2);
    }
    
    public static void updatePosition() 
    {
        //Position = Position + h*(Velocity)
        position.set(0, position.get(0) + (h * velocity.get(0)));
        position.set(1, position.get(1) + (h * velocity.get(1)));
        position.set(2, position.get(2) + (h * velocity.get(2)));   
    }
    
    public static void updateVelocity() {
        //Velocity = Velocity + h(acceleration)
        velocity.set(0, velocity.get(0) + (h * acceleration.get(0)));
        velocity.set(1, velocity.get(1) + (h * acceleration.get(1)));
        velocity.set(2, velocity.get(2) + (h * acceleration.get(2)));
    }
    
    public static void writeConditionsToFile() 
    {
        try(FileWriter fw = new FileWriter(output)) {
            PrintWriter print = new PrintWriter(new BufferedWriter(fw));
            if(output.length() == 0) {
                print.println("/----------- Initial Conditions -----------/");
                print.println("*************************************************\n");
            }
            print.print("Time: ");
            print.printf("%.2f", initialTime);
            print.print("s\n");
            print.println("Initial Normal: i = " + nBar.get(0) + ", j = " + nBar.get(1) + ", k = " + nBar.get(2) + ") m");
            print.println("Starting Position: (i = " + position.get(0) + ", j = " + position.get(1) + ", k = " + position.get(2) + ") m");
            print.println("Starting Velocity: (i = " + velocity.get(0) + ", j = " + velocity.get(1) + ", k = " + velocity.get(2) + ") m/s");
            print.println("Mu Kinetic: " + muKinetic);
            print.println("Mass: " + mass + "Kg");
            print.println("Area: " + area + "m^2");
            print.println("*************************************************\n\n");
            print.flush();
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    
    public static void writeStepToFile() {
        try(FileWriter fw = new FileWriter(output, true)) {
            PrintWriter print = new PrintWriter(new BufferedWriter(fw));
            if(isFirst) {
                print.println("/----------- Results -----------/");
                print.println("*************************************************\n");
                isFirst = false;
            }
            print.print("Time: ");
            print.printf("%.2f", currentTime);
            print.print("s\n");
            print.println("Current Mass: " + mass + "Kg");
            print.println("Current Force Gravity: i = " + forceGravity.get(0) + ", j = " + forceGravity.get(1) + ", k = " + forceGravity.get(2) + ") N");
            print.println("Current Force GN: i = " + forceGn.get(0) + ", j = " + forceGn.get(1) + ", k = " + forceGn.get(2) + ") N");
            print.println("Current Force GP: i = " + forceGp.get(0) + ", j = " + forceGp.get(1) + ", k = " + forceGp.get(2) + ") N");
            print.println("Current Frictional Force: i = " + forceFriction.get(0) + ", j = " + forceFriction.get(1) + ", k = " + forceFriction.get(2) + ") N");
            print.println("Current Force Normal: i = " + forceNormal.get(0) + ", j = " + forceNormal.get(1) + ", k = " + forceNormal.get(2) + ") m");
            print.println("Current Net Force: i = " + forceNet.get(0) + ", j = " + forceNet.get(1) + ", k = " + forceNet.get(2) + ") N");
            print.println("Current Acceleration: i = " + acceleration.get(0) + ", j = " + acceleration.get(1) + ", k = " + acceleration.get(2) + ") m/s^2");
            print.println("Position: i = " + position.get(0) + ", j = " + position.get(1) + ", k = " + position.get(2) + ") m");
            print.println("Velocity: i = " + velocity.get(0) + ", j = " + velocity.get(1) + ", k = " + velocity.get(2) + ") m/s");
            print.println("*************************************************\n");
            print.flush();
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    
    public static void initialize() {
        //More or less the same as last CA
        //NEEDS: time, radius, density, mass, position, velocity, normal force.
        System.out.print("Please enter the start time of the simulation in seconds: ");
        initialTime = in.nextDouble();
        System.out.print("Please enter the duration of the simulation in seconds: ");
        double duration = in.nextDouble();
        finalTime = initialTime + duration;
        
        System.out.println("Each step of the simulation will increment by 0.1 of a second.");
        
        System.out.print("Please enter the radius of sphere: ");
        double radius = in.nextDouble();
        System.out.print("Please enter the density of the sphere: ");
        double ballDensity = in.nextDouble();
        volume = 4.0/3 * Math.PI * (radius * radius * radius);
        mass = ballDensity * volume;
        area = Math.PI * (radius * radius);
        //Finds mass of the ball based on inputs for this specific simulation.
        
        System.out.println("");
        
        System.out.println("Please enter the three initial values for the normal force vector: ");
        System.out.print("i: " );
        double n1 = in.nextDouble();
        System.out.print("j: " );
        double n2 = in.nextDouble();
        System.out.print("k: " );
        double n3 = in.nextDouble();
        nBar.add(0, n1);
        nBar.add(1, n2);
        nBar.add(2, n3);
        System.out.println("");
        
        System.out.println("Please enter the three initial values for the position vector: ");
        System.out.print("i: " );
        double p1 = in.nextDouble();
        System.out.print("j: " );
        double p2 = in.nextDouble();
        System.out.print("k: " );
        double p3 = in.nextDouble();
        position.add(0, p1);
        position.add(1, p2);
        position.add(2, p3);
        System.out.println("");

        System.out.println("Please enter the three initial values for the velocity vector: ");
        System.out.print("i: " );
        double v1 = in.nextDouble();
        System.out.print("j: " );
        double v2 = in.nextDouble();
        System.out.print("k: " );
        double v3 = in.nextDouble();
        velocity.add(0, v1);
        velocity.add(1, v2);
        velocity.add(2, v3);
        System.out.println("");
        
        System.out.println("Please enter the values for muStatic and then muKinetic: ");
        System.out.print("muStatic: ");
        muStatic = in.nextDouble();
        System.out.print("muKinetic: ");
        muKinetic = in.nextDouble();
        
        writeConditionsToFile();
    }
}
