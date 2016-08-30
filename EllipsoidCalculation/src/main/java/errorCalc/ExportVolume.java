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
import biz.personalAcademics.ellipsoid.customExceptions.InvalidUserInputException;
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
			
			double volumeEst = ellip.getEstimatedVolume();
			System.out.println(i + "\n" + "Estimated: " + volumeEst);
			double volumeReal = ellip.getEstimatedVolume(10000000);
			System.out.println("Actual: " + volumeReal);
			System.out.println("Random: " + n);
			System.out.println("-----------------");
			double error = Math.abs(volumeEst - volumeReal);
			double errorUnit = error/volumeEst;
			
			double lowestEccentricity = 0;	
			double[] sortedAxis = ellip.getSortedAxes();
			double cOverA = ellip.getC()/ellip.getA();
			double bOverA = ellip.getB()/ellip.getA();
			
			if(cOverA < bOverA){
				lowestEccentricity = cOverA;
			}else{
				lowestEccentricity = bOverA;
			}
			
			fileOut.printf("%f,%f,%f\n", errorUnit, lowestEccentricity, ellip.getTotalEccentricity());
						
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
		
		public double getEstimatedVolume()
				throws InvalidUserInputException {
			if (sampleSize < Ellipsoid.MIN_SAMPLE_SIZE) {
				throw new InvalidUserInputException(
						new Integer(sampleSize).toString());
			}


				if (this.eccentricity >= 1) {
										
					return new EllipsoidSphericalCoords(startRadianTheta,
							endRadianTheta, radianMeasureOffZAxisStart,
							radianMeasureOffZAxisEnd, a, b, c)
							.getEstimatedVolume(sampleSize);

				} else {
					
					return new EllipsoidRectangularCoords(startRadianTheta,
							endRadianTheta, radianMeasureOffZAxisStart,
							radianMeasureOffZAxisEnd, a, b, c)
							.getEstimatedVolume(sampleSize);

				}

			
		}
	}

}
