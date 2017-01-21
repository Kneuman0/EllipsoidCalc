package errorCalc;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Random;

import static java.lang.Math.PI;
import biz.personalAcademics.ellipsoid.Ellipsoid;
import biz.personalAcademics.ellipsoid.EllipsoidRectangularCoords;
import biz.personalAcademics.ellipsoid.EllipsoidSphericalCoords;
import biz.personalAcademics.ellipsoidCalc.EllipsoidCalcController;
import biz.personalAcademics.lib.pathClasses.PathGetter;

public class ExportVolume {

	public static void main(String[] args) {
		PrintWriter fileOut = null;
		FileWriter file = null;
		try {
			String path = new PathGetter(new EllipsoidCalcController()).getAbsoluteSubfolderPath();
			file = new FileWriter(path + "/error.csv", true);
			fileOut = new PrintWriter(file);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		Random rand = new Random();
		
		for(int i = 0; i < 51; i++){
			double n = rand.nextDouble();
			double a = (rand.nextInt(300) + 1) * n;
			double b = (rand.nextInt(300) + 1) * n;
			double c = (rand.nextInt(300) + 1) * n;
			double[] sortedAxes = {a,b,c};
			Arrays.sort(sortedAxes);
			MonteCarloEllip ellip =  new MonteCarloEllip(0, 2*PI, 0, PI,
					sortedAxes[2], b, c);
			
			int sampleSize = 
//					rand.nextInt(50_000_000) + 1_000
					10_000_000;
			double volumeEst = ellip.getEstimatedVolume(sampleSize);
			System.out.println(i + "\n" + "Estimated: " + volumeEst);
			double volumeReal = ellip.getEstimatedVolume();
			System.out.println("Actual: " + volumeReal);
			System.out.println("Random: " + sampleSize);
			System.out.println("-----------------");
			double error = Math.abs(volumeEst - volumeReal);
			double errorUnit = error/volumeEst;
			
			double lowestEccentricity = 0;	
			double cOverA = ellip.getC()/ellip.getA();
			double bOverA = ellip.getB()/ellip.getA();
			
			if(cOverA < bOverA){
				lowestEccentricity = cOverA;
			}else{
				lowestEccentricity = bOverA;
			}
			
			fileOut.printf("errorUnit:,%f,SampleSize:,%d,lowestEccentricityTrace:,%.f \n",
					errorUnit, sampleSize, lowestEccentricity);
						
		}
		
		fileOut.close();
		try {
			file.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	public static class MonteCarloEllip extends Ellipsoid{
		
		public MonteCarloEllip(double startRadianTheta, double endRadianTheta,
				double radianMeasureOffZAxisStart, double radianMeasureOffZAxisEnd,
				double a, double b, double c) {
			super(startRadianTheta, endRadianTheta, radianMeasureOffZAxisStart,
				radianMeasureOffZAxisEnd, a, b, c);
		}
		
		public double getEstimatedVolume(int sampleSize){

				if (this.eccentricity >= 1) {
					
					this.sampleSize = sampleSize;
					
					// average 5 integrations
					double averageVolume = 0;
					for(int i = 0; i < 5; i++){
						averageVolume += new EllipsoidSphericalCoords(startRadianTheta,
								endRadianTheta, radianMeasureOffZAxisStart,
								radianMeasureOffZAxisEnd, a, b, c)
								.getEstimatedVolume(sampleSize);
					}
					
					averageVolume /= 5.0;
					
					return averageVolume;

				} else {

					this.sampleSize = sampleSize;
					
					double averageVolume = 0;
					for(int i = 0; i < 5; i++){
						averageVolume += new EllipsoidRectangularCoords(startRadianTheta,
								endRadianTheta, radianMeasureOffZAxisStart,
								radianMeasureOffZAxisEnd, a, b, c)
								.getEstimatedVolume(sampleSize);
					}
					
					averageVolume /= 5.0;
					
					return averageVolume;
				}
		}

	}
}
