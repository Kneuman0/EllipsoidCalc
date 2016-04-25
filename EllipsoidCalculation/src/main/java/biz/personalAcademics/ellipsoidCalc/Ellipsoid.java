package biz.personalAcademics.ellipsoidCalc;

import java.util.Random;


public class Ellipsoid {

	private double startRadianTheta, endRadianTheta, radianMeasureOffZAxisEnd,
			radianMeasureOffZAxisStart, a, b, c;
	
	Random portionOfAxis;
	

	public Ellipsoid(double startRadianTheta, double endRadianTheta,
			double radianMeasureOffZAxisEnd, double radianMeasureOffZAxisStart,
			double a, double b, double c) {

		this.startRadianTheta = startRadianTheta;
		this.endRadianTheta = endRadianTheta;
		this.radianMeasureOffZAxisEnd = radianMeasureOffZAxisEnd;
		this.radianMeasureOffZAxisStart = radianMeasureOffZAxisStart;
		this.a = a;
		this.b = b;
		this.c = c;
		
		/*
		 *  the nextDouble method will be called from this object and 
		 *  multiplied by one of the axes. This will yield a coordinate that lays 
		 *  within the 1st octant.
		 */
		
		portionOfAxis = new Random();
		
	}

	/**
	 * Finds the volume of the ellipsoid based off the spherical coordinates the
	 * user passed in
	 * 
	 * @param startRadianTheta
	 * @param endRadianTheta
	 * @param radianMeasureOffZAxisEnd
	 * @param radianMeasureOffZAxisStart
	 * @param a
	 * @param b
	 * @param c
	 * @return
	 */
	public double getExactVolume() {

		double oneThird = 1 / (double) 3;
		double volume = (-oneThird)
				* (a * b * c * (Math.cos(radianMeasureOffZAxisEnd) - Math
						.cos(radianMeasureOffZAxisStart)))
				* (endRadianTheta - startRadianTheta);
		return volume;
	}
	
	public double getEstimatedVolume(int sampleSize){
		int insideShape = 0;
		
		for(int i = 0; i < 500; i++){
			if(determineIfPointIsInsideShape(generate3DSamplePoint(), startRadianTheta, 
					endRadianTheta, radianMeasureOffZAxisStart, radianMeasureOffZAxisEnd)){
				insideShape++;
			}
			
			
		}
		
		double portionOfKnownVolume = insideShape / (double)sampleSize;
		
		return portionOfKnownVolume * a * b * c;
	}
	
//	private boolean determineXCoordinatePosition(double[] coord){
//		double xLimit = a * Math.sqrt(1 - Math.pow(coord[2], 2)/(c * c) 
//				- Math.pow(coord[1], 2)/(b * b));
//		
//		/*
//		 *  If the generated point is less than the curve of the function then 
//		 *  the point lays within the volume of the ellipsoid, otherwise its isn't.
//		 */
//		return coord[0] <= xLimit;
//	}
//	
//	private boolean determineYCoordinatePosition(double[] coord){
//		double yLimit = b * Math.sqrt(1 - Math.pow(coord[2], 2)/(c * c) 
//				- Math.pow(coord[0], 2)/(a * a));
//		
//		/*
//		 *  If the generated point is less than the curve of the function then 
//		 *  the point lays within the volume of the ellipsoid, otherwise its isn't.
//		 */
//		return coord[1] <= yLimit;
//	}
	
	private boolean determineIfPointIsInsideShape(double[] coord, double thetaStart, double thetaEnd,
			double phiStart, double phiEnd){
		boolean insideShape = false;
		
		determineXAxisToThetaPhiLineToXYPlane(coord, thetaStart, thetaEnd, phiStart, phiEnd);
		
		/*
		 *  MUST DETERMINE ALL OTHER SCENERIOS AND CODE SOLUTIONS TO EACH!
		 */
		return insideShape;
	}
	
	
	/**
	 * Determines if the point is between the x axis and the line y = tan(theta) * x
	 * Determines if the point is between the line z = tan(phi) * y
	 * Determines if the point is inside the equation of the ellipsoid.
	 * @param coord
	 * @param thetaStart
	 * @param thetaEnd
	 * @param phiStart
	 * @param phiEnd
	 * @return
	 */
	private boolean determineXAxisToThetaPhiLineToXYPlane(double[] coord, double thetaStart, double thetaEnd,
			double phiStart, double phiEnd){
		
		boolean insideShape = false;
		
		if(determineIfPointIsBetweenPhiStartAndXYPlane(coord, phiStart) ||
			determineIfPointIsBetweenXAxisAndThetaStart(coord, thetaEnd) ||
					determineIfInsideEllipsoid(coord)){
			
			insideShape = true;
				
			}
		return insideShape;
	}
	
	
	/**
	 * This method determines whether or not the point lays between thetaStart
	 * and thetaEnd. This method is ONLY used when thetaEnd is NOT pi/2
	 * @param coord
	 * @param theta
	 * @return
	 */
	private boolean determineIfBetweenThetaValues(double[] coord, double thetaStart, double thetaEnd){
		boolean insideShape = false;
		
		double yLower = Math.tan(thetaStart) * coord[0];
		double yUpper = Math.tan(thetaEnd) * coord[0];
		
		if(coord[1] <= yUpper || yLower <= coord[1]){
			insideShape = true;
		}
		
		return insideShape;
	}
	
	/**
	 * Determines whether or not point is between the line y = tan(theta) * x and the y axis.
	 * @param coord
	 * @param thetaStart
	 * @return
	 */
	private boolean determineIfPointIsBetweenYAxisAndThetaStart(double[] coord, double thetaStart){
		boolean insideShape = false;
		
		double lowerYLimit = Math.tan(thetaStart) * coord[0];
		
		if(lowerYLimit <= coord[1]){
			insideShape = true;
		}
		
		return insideShape;
	}
	
	/**
	 * Determines whether or not point is between the line y = tan(theta) * x and the x axis.
	 * @param coord
	 * @param thetaStart
	 * @return
	 */
	private boolean determineIfPointIsBetweenXAxisAndThetaStart(double[] coord, double thetaEnd){
		boolean insideShape = false;
		
		double upperYLimit = Math.tan(thetaEnd) * coord[0];
		
		if(coord[1] <= upperYLimit){
			insideShape = true;
		}
		
		return insideShape;
	}
	
	
	/**
	 * Use this method ONLY if phiStart is not 0 and phiEnd is NOT pi/2
	 * @param coord
	 * @param phiStart
	 * @param phiEnd
	 * @return
	 */
	private boolean determineIfPointIsBetweenPhiCones(double[] coord, double phiStart, double phiEnd){
		boolean insideShape = false;
		
		if(coord[2] <= getZValueOfCone(coord, phiStart) || getZValueOfCone(coord, phiEnd) <= coord[2]){
			insideShape = true;
		}
		
		return insideShape;
	}
	
	private boolean determineIfPointIsBetweenPhiStartAndXYPlane(double[] coord, double phiStart){
		boolean insideShape = false;
		
		if(coord[2] <= getZValueOfCone(coord, phiStart)){
			insideShape = true;
		}
		
		return insideShape;
		
	}
	
	private boolean determineIfPointIsBetweenXZPlaneAndCone(double[] coord, double phiEnd){
		boolean insideShape = false;
		
		if(getZValueOfCone(coord, phiEnd) <= coord[2]){
			insideShape = true;
		}
		
		return insideShape;
	}
	
	private boolean determineIfInsideEllipsoid(double[] coord){
		boolean insideShape = false;
		
		double zLimit = c * Math.sqrt(1 - Math.pow(coord[0], 2)/(a * a) 
				- Math.pow(coord[1], 2)/(b * b));
		/*
		 * Checks if the point is inside the function. If all previous 
		 * checks are true, and this check is true, the point is inside
		 * of the shape.
		 */
		if(coord[2] <= zLimit){
			insideShape = true;
		}
		
		return insideShape;
	}
	
	private double getZValueOfCone(double[] coord, double phi){
		return Math.sqrt((coord[0] * coord[0]) + (coord[1] * coord[1]))/ Math.tan(phi);
	}
	
	/**
	 * generates a 3 dimensional sample point that lays with the first
	 * octant of the ellipsoid Cartesian field.
	 * @return
	 */
	public double[] generate3DSamplePoint(){
		double[] samplePoint = new double[3];
		
		samplePoint[0] =  a * portionOfAxis.nextDouble(); 
		samplePoint[1] = b * portionOfAxis.nextDouble();
		samplePoint[2] = c * portionOfAxis.nextDouble();
		
		return samplePoint;
	}

	/**
	 * @return the startRadianTheta
	 */
	public double getStartRadianTheta() {
		return startRadianTheta;
	}

	/**
	 * @param startRadianTheta
	 *            the startRadianTheta to set
	 */
	public void setStartRadianTheta(double startRadianTheta) {
		this.startRadianTheta = startRadianTheta;
	}

	/**
	 * @return the endRadianTheta
	 */
	public double getEndRadianTheta() {
		return endRadianTheta;
	}

	/**
	 * @param endRadianTheta
	 *            the endRadianTheta to set
	 */
	public void setEndRadianTheta(double endRadianTheta) {
		this.endRadianTheta = endRadianTheta;
	}

	/**
	 * @return the radianMeasureOffZAxisEnd
	 */
	public double getRadianMeasureOffZAxisEnd() {
		return radianMeasureOffZAxisEnd;
	}

	/**
	 * @param radianMeasureOffZAxisEnd
	 *            the radianMeasureOffZAxisEnd to set
	 */
	public void setRadianMeasureOffZAxisEnd(double radianMeasureOffZAxisEnd) {
		this.radianMeasureOffZAxisEnd = radianMeasureOffZAxisEnd;
	}

	/**
	 * @return the radianMeasureOffZAxisStart
	 */
	public double getRadianMeasureOffZAxisStart() {
		return radianMeasureOffZAxisStart;
	}

	/**
	 * @param radianMeasureOffZAxisStart
	 *            the radianMeasureOffZAxisStart to set
	 */
	public void setRadianMeasureOffZAxisStart(double radianMeasureOffZAxisStart) {
		this.radianMeasureOffZAxisStart = radianMeasureOffZAxisStart;
	}

	/**
	 * @return the a
	 */
	public double getA() {
		return a;
	}

	/**
	 * @param a
	 *            the a to set
	 */
	public void setA(double a) {
		this.a = a;
	}

	/**
	 * @return the b
	 */
	public double getB() {
		return b;
	}

	/**
	 * @param b
	 *            the b to set
	 */
	public void setB(double b) {
		this.b = b;
	}

	/**
	 * @return the c
	 */
	public double getC() {
		return c;
	}

	/**
	 * @param c
	 *            the c to set
	 */
	public void setC(double c) {
		this.c = c;
	}

	/**
	 * Converts value in text to a decimal
	 * 
	 * @param input
	 * @return
	 * @throws InvalidUserInputException
	 */
	public static double convertToDecimal(String input)
			throws InvalidUserInputException {
		String[] fraction = input.split("/");
		// Returns user input if array is not a fraction
		if (fraction.length == 1) {
			try {
				return Double.parseDouble(input);
			} catch (NumberFormatException e) {
				throw new InvalidUserInputException(input);
			}
			// converts user input into a decimal if a fraction or throws
			// invalid input exception
		} else {
			int[] numbers = new int[2];
			for (int i = 0; i < fraction.length; i++) {
				try {
					numbers[i] = Integer.parseInt(fraction[i]);
				} catch (NumberFormatException | ArrayIndexOutOfBoundsException e) {
					throw new InvalidUserInputException(input);
				}
			}
			return numbers[0] / (double) numbers[1];
		}

	}

	/**
	 * If user inputs angle in degrees, it is converted and returned in radians.
	 * If radians were chosen, pi is multiplied in
	 * 
	 * @param degrees
	 * @param degreeIsSelected
	 * @return
	 * @throws InvalidUserInputException
	 */
	public static double convertThetaToRadians(String degrees,
			boolean degreeIsSelected) throws InvalidUserInputException {
		double testDegree;

		try {
			testDegree = Double.parseDouble(degrees);
		} catch (NumberFormatException e1) {
			throw new InvalidUserInputException(degrees);
		}

		if (degreeIsSelected) {
			return (testDegree / 180) * (Math.PI);
		} else {
			return testDegree * Math.PI;
		}
	}

	public String toString() {
		return String.format("%f", this.getExactVolume());
	}

}
