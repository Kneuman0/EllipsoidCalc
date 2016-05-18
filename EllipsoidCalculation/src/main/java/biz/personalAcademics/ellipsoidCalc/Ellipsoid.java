package biz.personalAcademics.ellipsoidCalc;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

public class Ellipsoid {

	private double startRadianTheta, endRadianTheta, radianMeasureOffZAxisEnd,
			radianMeasureOffZAxisStart, a, b, c;

	private final int x = 0, y = 1, z = 2;
	private final int phi = 0, theta = 1, p = 2;
	private final int r = 0;
	
	private double[] sortedAxes;

	private double sampleSize;

	private Random randomGenerator;

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
		
		sortedAxes = new double[]{a, b, c};
		Arrays.sort(sortedAxes);

		/*
		 * the nextDouble method will be called from this object and multiplied
		 * by one of the axes. This will yield a coordinate that lays within the
		 * 1st octant.
		 */

		randomGenerator = new Random();

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

	/**
	 * Estimates the volume of the ellipsoid using the specified number of
	 * random sample points and the most reliable method developed thus far
	 * 
	 * @param sampleSize
	 * @return
	 */
	public double getEstimatedVolume(int sampleSize) {
		return getEstimatedVolumeRect(sampleSize);
	}

	/**
	 * Estimates the volume of the portion of the ellipsoid using using the
	 * most reliable method developed thus far
	 * 
	 * @return
	 */
	public double getEstimatedVolume() {
		return getEstimatedVolumeRect();
	}
	
	public double getEstimatedVolumeCylinder(){
		int insideShape = 0;
		final int sampleSize = 5_000_000;

		for (int i = 0; i < sampleSize; i++) {
			if (pointInsideShapeCyl(generate3DSamplePointCyl())) {
				insideShape++;
			}

		}
		
		double portionOfKnownVolume = insideShape / (double) sampleSize;
		System.out.printf("\nportion = %.3f, longest Axis = %.3f", portionOfKnownVolume, sortedAxes[2]);
		
		// portion of volume of 3D system V = (pi)(h)(r^2)
		return portionOfKnownVolume * Math.PI * Math.pow(getLongestAxis(), 2) * 2 * this.c;
	}
	
	/**
	 * creates a specified number of sample points inside a sphere using a spherical coordinate system
	 * the sphere's volume is defined as (4/3)pi(largestAxis)^3. This method inscribes an
	 * ellipsoid inside the sphere to compute the estimated volume. The error has
	 * empirically been found to be as much as +/- .03
	 * @return
	 */
	public double getEstimatedVolumeSphere(){
		int insideShape = 0;
		final int sampleSize = 10_000_000;

		for (int i = 0; i < sampleSize; i++) {
			if (determineIfPointIsInsideShapeRandomDistributionSphere(generateRandom3DSamplePointSphere())) {
				insideShape++;
			}

		}
		
		double portionOfKnownVolume = insideShape / (double) sampleSize;
		System.out.printf("\nportion = %.3f, longest Axis = %.3f", portionOfKnownVolume, sortedAxes[2]);
		
		// portion of volume of 3D system
		return portionOfKnownVolume * (4.0/3.0) * Math.PI * Math.pow(sortedAxes[2], 3);
	}
	
	/**
	 * creates 30 million sample points inside a sphere using a spherical coordinate system
	 * the sphere's volume is defined as (4/3)pi(largestAxis)^3. This method inscribes an
	 * ellipsoid inside the sphere to compute the estimated volume. The error has
	 * empirically been found to be as much as +/- .03
	 * @return
	 */
	public double getEstimatedVolumeSphere(int sampleSize){
		int insideShape = 0;

		for (int i = 0; i < sampleSize; i++) {
			if (determineIfPointIsInsideShapeRandomDistributionSphere(generateRandom3DSamplePointSphere())) {
				insideShape++;
			}

		}
		
		double portionOfKnownVolume = insideShape / (double) sampleSize;
		
		// portion of volume of 3D system
		return portionOfKnownVolume * (4.0/3.0) * Math.PI * Math.pow(sortedAxes[2], 3);
	}
	
	/**
	 * creates 10 million sample points inside a rectangular prism using a rectangular coordinate system
	 * the rectangular prism's volume is defined as 2a * 2b * 2c where a b and c represent the 
	 * axes of the ellipsoid. This method inscribes an
	 * ellipsoid inside the rectangular prism to compute the estimated volume. The error has
	 * empirically been found to be as much as +/- .1
	 * @return
	 */
	public double getEstimatedVolumeRect(){
		int insideShape = 0;
		final int sampleSize = 10_000_000;

		for (int i = 0; i < sampleSize; i++) {
			if (determineIfPointIsInsideShapeRandomDistributionRect(generateRandom3DSamplePointRect())) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) sampleSize;

		// portion of volume of 3D system
		return portionOfKnownVolume * (2 * a) * (2 * b) * (2 * c);
	}
	
	/**
	 * creates a specified number of sample points inside a rectangular prism using a rectangular coordinate system
	 * the rectangular prism's volume is defined as 2a * 2b * 2c where a b and c represent the 
	 * axes of the ellipsoid. This method inscribes an
	 * ellipsoid inside the rectangular prism to compute the estimated volume. The error has
	 * empirically been found to be as much as +/- .1
	 * @return
	 */
	public double getEstimatedVolumeRect(int sampleSize){
		int insideShape = 0;

		for (int i = 0; i < sampleSize; i++) {
			if (determineIfPointIsInsideShapeRandomDistributionRect(generateRandom3DSamplePointRect())) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) sampleSize;

		// portion of volume of 3D system
		return portionOfKnownVolume * (2 * a) * (2 * b) * (2 * c);
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
	
	private double getLongestAxis(){
		return sortedAxes[2];
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

			// if in fraction form, it is converted to decimal, otherwise
			// invalidUserInputException is thrown
			testDegree = Ellipsoid
					.convertToDecimal(degrees.replaceAll(" ", ""));
		}

		if (degreeIsSelected) {
			return (testDegree / 180) * (Math.PI);
		} else {
			return testDegree * Math.PI;
		}
	}

	/**
	 * This method will return the estimated volume of the shape through random
	 * sample points rounded to two decimal places.
	 */
	public String toString() {
		return String.format("%.2f +/- .05", this.getEstimatedVolume());
	}
	
	//-------------------Random Dist. Cylinder------------------------------
	
	private boolean pointInsideShapeCyl(double[] coord){
		
		if(checkIfThetaIsBetweenStartAndEndCyl(coord) &&
				pointBetweenPhiValuesCyl(coord) &&
				pointInsideEquationOfEllipsoidCyl(coord) &&
				checkLengthOfR(coord)){
			
			return true;
			
		}else{
			
			return false;
		}
	}
	
	private boolean checkLengthOfR(double[] coord){
		double numerator = this.a * this.b;
		
		// ( b^2 cos^2(theta) + a^2 sin^2(theta) )^(1/2)
		double denominator = Math.sqrt(Math.pow(this.b * Math.cos(coord[theta]), 2)
				+ Math.pow(this.a * Math.sin(coord[theta]), 2) );
		double rValue = numerator / denominator;
		
		return coord[r] <= rValue;
		
	}
	
	private boolean pointInsideEquationOfEllipsoidCyl(double[] coord){
		double xValue = Math.pow((coord[r] * Math.cos(coord[theta])) / this.a, 2);
		double yValue = Math.pow((coord[r] * Math.sin(coord[theta])) / this.b, 2);
		double zValue = Math.pow(coord[z]/this.c, 2);
		
		double pointValue = xValue + yValue + zValue;
		
		if(pointValue <= 1){
			return true;
		}else{
			return false;
		}
	}
	
	private boolean pointBetweenPhiValuesCyl(double[] coord){
		
		if(checkThatZMoreThanPhiLineEnd(coord) && 
				checkThatZLessThanPhiLineStart(coord)){
			return true;
		}else{
			return false;
		}
	}
	
	private boolean checkThatZMoreThanPhiLineEnd(double[] coord){
		if(radianMeasureOffZAxisEnd == Math.PI){
			
			return true;
		}else{
			
			double phiValue = ( Math.cos(radianMeasureOffZAxisEnd) /
					Math.sin(radianMeasureOffZAxisEnd) ) * coord[r];
			
			if(phiValue <= coord[z]){
				
				return true;
				
			}else{
				
				return false;
			}
		}
	}
	
	
	private boolean checkThatZLessThanPhiLineStart(double[] coord){
		if(radianMeasureOffZAxisStart == 0){
			
			return true;
		}else{
			
			double phiValue = ( Math.cos(radianMeasureOffZAxisStart) /
					Math.sin(radianMeasureOffZAxisStart) ) * coord[r];
			
			if(coord[z] <= phiValue){
				
				return true;
				
			}else{
				
				return false;
			}
		}
	}
	
	private boolean checkIfThetaIsBetweenStartAndEndCyl(double[] coord){
		boolean insideThetaBound = false;

		if (startRadianTheta <= coord[theta]
				&& coord[theta] <= endRadianTheta) {

			insideThetaBound = true;
		}
	
		return insideThetaBound;
	}
	
	private double[] generate3DSamplePointCyl(){
		double[] samplePoint = new double[3];
		
		// r will range from 0 to the longest axis
		samplePoint[r] = getLongestAxis() * Math.sqrt(randomGenerator.nextDouble());
						
		// theta will range 0 to 2pi
		samplePoint[theta] = 2 * Math.PI * randomGenerator.nextDouble();
		
		// z will range from 0 to the c axis both above and below the xy plane
		samplePoint[z] = this.c * randomGenerator.nextDouble() * getRandomNegation();
		
		return samplePoint;		
	}
	
	//-------------------Random Dist. Sphere---------------------------------
	
	/**
	 * This method can only be used for Spherical coordinates. Checks if point is
	 * between specified bounds for phi. Checks if point is within the equation
	 * of an ellipsoid. 
	 * 
	 * fomulas: theta: thetaStart < thetaSample < thetaEnd phi: phiStart <
	 * phiSample < phiEnd equation of ellipsoid: (x^2/a^2) + (y^2/b^2) +
	 * (z^2/c^2) <= 1
	 * 
	 * @param coord
	 * @return
	 */
	private boolean determineIfPointIsInsideShapeRandomDistributionSphere(
			double[] coord) {

		boolean insideShape = false;

		if (pointInsideEquationOfEllipsoidSphere(coord)
				&& pointBetweenPhiStartAndPhiEndSphere(coord)
				&& pointBetweenThetaStartAndThetaEndRandomDistSphere(coord)) {
			insideShape = true;
		}

		return insideShape;
	}
	
	/**
	 * Checks if the phi value of sample point is >= phiStart but <= phiEnd
	 * 
	 * @param coord
	 * @return
	 */
	private boolean pointBetweenPhiStartAndPhiEndSphere(double[] coord) {
		boolean insidePhiBound = false;

			if (radianMeasureOffZAxisStart <= coord[phi]
					&& coord[phi] <= radianMeasureOffZAxisEnd) {

				insidePhiBound = true;
			}
		
		return insidePhiBound;
	}

	/**
	 * checks if theta value of sample point is > thetaStart but < thetaEnd.
	 * Randomly decides whether or not to include special cases.
	 * 
	 * @param coord
	 * @return
	 * @throws NegativeZAxisException
	 */
	private boolean pointBetweenThetaStartAndThetaEndRandomDistSphere(double[] coord) {
		boolean insideThetaBound = false;

			if (startRadianTheta <= coord[theta]
					&& coord[theta] <= endRadianTheta) {

				insideThetaBound = true;
			}
		

		return insideThetaBound;
	}
	
	/**
	 * checks if sample point is inside the equation of an ellipsoid using the
	 * boolean expression: 
	 * ((psin(phi)cos(theta))^2/a^2) + ((psin(phi)cos(theta))^2/b^2) + 
	 * ((pcos(phi))^2/c^2) <= 1
	 * 
	 * @param coord
	 * @return
	 */
	private boolean pointInsideEquationOfEllipsoidSphere(double[] coord) {
		boolean insideShape = false;
		
		double xSquared = Math.pow((coord[p] * Math.sin(coord[phi]) * Math.cos(coord[theta])), 2);
		double ySquared = Math.pow((coord[p] * Math.sin(coord[phi]) * Math.sin(coord[theta])), 2);
		double zSquared = Math.pow((coord[p] * Math.cos(coord[phi])), 2);
		
		/*
		
		point is inside shape if ((p * sin(phi)cos(theta))^2/a^2) + ((p * sin(phi)cos(theta))^2/b^2)
	 	+ ((p * cos(phi))^2/c^2) <= 1 then it is inside the ellipsoid
		 
		*/
		double valueOfPointInEllipsoidEquation = (xSquared)
				/ (a * a) + (ySquared) / (b * b)
				+ (zSquared) / (c * c);

		if (valueOfPointInEllipsoidEquation <= 1) {
			insideShape = true;
		}

		return insideShape;
	}
	
	/**
	 * generates a 3 dimensional sample point using spherical coordinates
	 * that lies somewhere inside the sphere (4/3)pi(r^2) where r is the 
	 * longest axis
	 * 
	 * @return
	 */
	private double[] generateRandom3DSamplePointSphere() {
		double[] samplePoint = new double[3];

		// phi ranges from 0 to pi
		samplePoint[phi] = Math.acos(2 * randomGenerator.nextDouble() - 1);
		
		// theta ranges from 0 to 2pi
		samplePoint[theta] = 2 * Math.PI * randomGenerator.nextDouble();
		
		// p ranges from 0 to the longest axis
		samplePoint[p] = sortedAxes[2] * Math.cbrt(randomGenerator.nextDouble());

		return samplePoint;
	}

	
	//--------------------Random Dist. Rectangular Prism----------------------
	
	/**
	 * This method can be used on uniformly distributed sample points. For
	 * random distribution use 'determineIfPointIsInsideShapeRandomDistribution'
	 * Checks if point is between specified bounds for theta. Checks if point is
	 * between specified bounds for phi. Checks if point is within the equation
	 * of an ellipsoid. Throws exceptions for cases where the point is on the
	 * edge of the specified shape.
	 * 
	 * fomulas: theta: thetaStart < thetaSample < thetaEnd phi: phiStart <
	 * phiSample < phiEnd equation of ellipsoid: (x^2/a^2) + (y^2/b^2) +
	 * (z^2/c^2) <= 1
	 * 
	 * @param coord
	 * @return
	 */
	private boolean determineIfPointIsInsideShapeRandomDistributionRect(
			double[] coord) {

		boolean insideShape = false;

		if (pointInsideEquationOfEllipsoidRect(coord)
				&& pointBetweenPhiStartAndPhiEndRect(coord)
				&& pointBetweenThetaStartAndThetaEndRandomDistRect(coord)) {
			insideShape = true;
		}

		return insideShape;
	}

	/**
	 * checks if sample point is inside the equation of an ellipsoid using the
	 * boolean expression: (x^2/a^2) + (y^2/b^2) + (z^2/c^2) <= 1
	 * 
	 * @param coord
	 * @return
	 */
	private boolean pointInsideEquationOfEllipsoidRect(double[] coord) {
		boolean insideShape = false;

		// point is inside shape if x^2/a^2 + y^2/b^2 + z^2/c^2 <= 1 then it is
		// inside the ellipsoid
		double valueOfPointInEllipsoidEquation = (coord[x] * coord[x])
				/ (a * a) + (coord[y] * coord[y]) / (b * b)
				+ (coord[z] * coord[z]) / (c * c);

		if (valueOfPointInEllipsoidEquation <= 1) {
			insideShape = true;
		}

		return insideShape;
	}

	/**
	 * Checks if the phi value of sample point is >= phiStart but <= phiEnd
	 * 
	 * @param coord
	 * @return
	 */
	private boolean pointBetweenPhiStartAndPhiEndRect(double[] coord)
			throws PointOnEdgeException {
		boolean insidePhiBound = false;

		try {
			if (radianMeasureOffZAxisStart <= getPhiValueOfSampleCoord(coord)
					&& getPhiValueOfSampleCoord(coord) <= radianMeasureOffZAxisEnd) {

				if (radianMeasureOffZAxisStart == getPhiValueOfSampleCoord(coord)
						|| getPhiValueOfSampleCoord(coord) == radianMeasureOffZAxisEnd) {

					throw new PointOnEdgeException(coord);
				}

				insidePhiBound = true;

			}
		} catch (PointOnOriginException e) {
			insidePhiBound = true;
			System.out.println(e.getMessage());
		}
		return insidePhiBound;
	}

	/**
	 * checks if theta value of sample point is > thetaStart but < thetaEnd.
	 * Randomly decides whether or not to include special cases.
	 * 
	 * @param coord
	 * @return
	 * @throws NegativeZAxisException
	 */
	private boolean pointBetweenThetaStartAndThetaEndRandomDistRect(double[] coord) {
		boolean insideThetaBound = false;

		try {
			if (startRadianTheta <= getThetaValueOfSampleCoordStartTheta(coord)
					&& getThetaValueOfSampleCoordEndTheta(coord) <= endRadianTheta) {

				insideThetaBound = true;
			}
		} catch (NegativeArraySizeException | PointOnOriginException
				| PositiveZAxisException e) {
			System.out.println(e.getMessage());
			// Randomly determines whether or not to include point in shape
			insideThetaBound = randomGenerator.nextBoolean();
		}

		return insideThetaBound;
	}

	/**
	 * Calculates the value of phi using the formula: arcCos(z/sqrt(x^2 + y^2 +
	 * z^2)). If z is negative, pi/2 is added to the angle. If z = 0 but x != 0
	 * or y != 0 then phi is pi/2. If x == 0 and y == 0 and z != 0, then phi is
	 * either 0 or pi. If x == 0, y == 0, z == 0, then PointOnOrigin exception
	 * is thrown.
	 * 
	 * @param coord
	 * @return
	 * @throws PointOnOriginException
	 */
	private double getPhiValueOfSampleCoord(double[] coord)
			throws PointOnOriginException {
		// above xy plane (positive z value)
		if (coord[z] > 0 && coord[x] != 0 && coord[y] != 0) {

			// arcCos(z/sqrt(x^2 + y^2 + z^2))
			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])));

			// below xy plane (negative z value)
		} else if (coord[z] < 0 && coord[x] != 0 && coord[y] != 0) {

			// arcCos(z/sqrt(x^2 + y^2 + z^2)) + pi
			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])))
			// + Math.PI
			;

			// on positive z axis = 0
		} else if (coord[z] > 0 && coord[x] == 0 && coord[y] == 0) {
			return 0;

			// on negative z axis = pi
		} else if (coord[z] < 0 && coord[x] == 0 && coord[y] == 0) {
			return Math.PI;

			// in zy plane above xy plane
		} else if (coord[z] > 0 && coord[x] == 0 && coord[y] != 0) {

			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])));

			// in zy plane below xy plane
		} else if (coord[z] < 0 && coord[x] == 0 && coord[y] != 0) {

			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])))
			// + Math.PI
			;

			// in zx plane above xy plane
		} else if (coord[z] > 0 && coord[x] != 0 && coord[y] == 0) {

			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])));

			// in zx plane below xy plane
		} else if (coord[z] < 0 && coord[x] != 0 && coord[y] == 0) {

			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])));

			// in xy plane
		} else if (coord[z] == 0) {

			if (coord[x] != 0 || coord[y] != 0) {
				return Math.PI / 2;

			} else {

				throw new PointOnOriginException(coord);

			}

			// on origin
		} else {

			throw new PointOnOriginException(coord);

		}

		// do case where phi is below x y plane

	}

	/**
	 * Returns the values of theta in the x y plane using the formula
	 * arcTan(y/x). ***If point is on positive x axis then 0 is returned.*** If
	 * point is in first quadrant then the angle is returned. If point is in
	 * second quadrant, then pi/2 is added to angle. If point is in third
	 * quadrant, then pi/2 is added to angle. If point is in fourth quadrant
	 * then 2pi is added to the angle. If x = 0, y = 0 and z > 0, then
	 * PositiveZAxisException is thrown. If x = 0, y = 0 and z < 0 then
	 * NegativeArraySizeException is thrown. If x = 0, y = 0 and z = 0, then
	 * PointOnOriginException is thrown.
	 * 
	 * @param coord
	 * @return
	 * @throws PositiveZAxisException
	 * @throws NegativeArraySizeException
	 * @throws PointOnOriginException
	 */
	private double getThetaValueOfSampleCoordStartTheta(double[] coord)
			throws PositiveZAxisException, NegativeArraySizeException,
			PointOnOriginException {
		// first quadrant
		if (coord[x] > 0 && coord[y] > 0) {
			return Math.atan(coord[y] / coord[x]);

			// Second quadrant
		} else if (coord[x] < 0 && coord[y] > 0) {
			return Math.atan(coord[y] / coord[x]) + Math.PI;

			// Third Quadrant
		} else if (coord[x] < 0 && coord[y] < 0) {
			return Math.atan(coord[y] / coord[x]) + Math.PI;

			// Fourth Quadrant
		} else if (coord[x] > 0 && coord[y] < 0) {
			return Math.atan(coord[y] / coord[x]) + 2 * Math.PI;

			// on positive y axis
		} else if (coord[x] == 0 && coord[y] > 0) {
			return Math.PI / 2;

			// on negative y axis
		} else if (coord[x] == 0 && coord[y] < 0) {
			return (3 / 2) * Math.PI;

			// on positive x axis (choosing 0 for radian measure because this is
			// the start value)
		} else if (coord[x] > 0 && coord[y] == 0) {
			return 0;

			// on negative x axis
		} else if (coord[x] < 0 && coord[y] == 0) {
			return Math.PI;

			// on negative z axis
		} else if (coord[x] == 0 && coord[y] == 0 && coord[z] < 0) {

			throw new NegativeZAxisException(coord);

		} else if (coord[x] == 0 && coord[y] == 0 && coord[z] > 0) {
			throw new PositiveZAxisException(coord);

			// on positive z axis
		} else {
			throw new PointOnOriginException(coord);
		}
	}

	/**
	 * Returns the values of theta in the x y plane using the formula
	 * arcTan(y/x). ***If point is on positive x axis then 2pi is returned.***
	 * If point is in first quadrant then the angle is returned. If point is in
	 * second quadrant, then pi/2 is added to angle. If point is in third
	 * quadrant, then pi/2 is added to angle. If point is in fourth quadrant
	 * then 2pi is added to the angle. If x = 0, y = 0 and z > 0, then
	 * PositiveZAxisException is thrown. If x = 0, y = 0 and z < 0 then
	 * NegativeArraySizeException is thrown. If x = 0, y = 0 and z = 0, then
	 * PointOnOriginException is thrown.
	 * 
	 * @param coord
	 * @return
	 * @throws PositiveZAxisException
	 * @throws NegativeArraySizeException
	 * @throws PointOnOriginException
	 */
	private double getThetaValueOfSampleCoordEndTheta(double[] coord)
			throws PositiveZAxisException, NegativeArraySizeException,
			PointOnOriginException {
		// first quadrant
		if (coord[x] > 0 && coord[y] > 0) {
			return Math.atan(coord[y] / coord[x]);

			// Second quadrant
		} else if (coord[x] < 0 && coord[y] > 0) {
			return Math.atan(coord[y] / coord[x]) + Math.PI;

			// Third Quadrant
		} else if (coord[x] < 0 && coord[y] < 0) {
			return Math.atan(coord[y] / coord[x]) + Math.PI;

			// Fourth Quadrant
		} else if (coord[x] > 0 && coord[y] < 0) {
			return Math.atan(coord[y] / coord[x]) + 2 * Math.PI;

			// on positive y axis
		} else if (coord[x] == 0 && coord[y] > 0) {
			return Math.PI / 2;

			// on negative y axis
		} else if (coord[x] == 0 && coord[y] < 0) {
			return (3 / 2) * Math.PI;

			// on positive x axis (choosing 0 for radian measure because this is
			// the start value)
		} else if (coord[x] > 0 && coord[y] == 0) {
			return 2 * Math.PI;

			// on negative x axis
		} else if (coord[x] < 0 && coord[y] == 0) {
			return Math.PI;

			// on negative z axis
		} else if (coord[x] < 0 && coord[y] == 0 && coord[z] < 0) {

			throw new NegativeZAxisException(coord);

			// on positive z axis
		} else if (coord[x] < 0 && coord[y] == 0 && coord[z] > 0) {

			throw new PositiveZAxisException(coord);

			// on origin
		} else {
			throw new PointOnOriginException(coord);
		}
	}

	/**
	 * Returns a positive 1 or -1 randomly
	 * 
	 * @return
	 */
	private int getRandomNegation() {
		boolean negativeSign = randomGenerator.nextBoolean();

		if (negativeSign) {
			return 1;
		} else {
			return -1;
		}
	}

	/**
	 * generates a 3 dimensional sample point that lays somewhere within the
	 * three dimensional coordinate system where each axis ranges from x = +/-
	 * a, y = +/- b and z = +/- c
	 * 
	 * @return
	 */
	private double[] generateRandom3DSamplePointRect() {
		double[] samplePoint = new double[3];

		samplePoint[x] = a * randomGenerator.nextDouble() * getRandomNegation();
		samplePoint[y] = b * randomGenerator.nextDouble() * getRandomNegation();
		samplePoint[z] = c * randomGenerator.nextDouble() * getRandomNegation();

		return samplePoint;
	}


	// --------------------Uniform distribution Monte Carlo--------------------
	/**
	 * Uses a uniformly distributed set of 3 dimentional coordinates to estimate
	 * volume. Formula for number of points is = [(number of sample points on 1
	 * positive axis)^3] * 8. Due to computation time, samples are limited to
	 * ~2,700,000 sample points.
	 * 
	 * @return
	 */
	public double getEstimatedVolumeThruUniformDistributionMonteCarlo() {
		double insideShape = 0;
		ArrayList<double[]> coord = generateUniformDistributionOfSamplePoints();

		for (int i = 0; i < coord.size(); i++) {
			try {

				if (determineIfPointIsInsideShapeNonRandomDistribution(coord
						.get(i))) {
					insideShape += 1;
				}

			} catch (PointOnEdgeException | NegativeZAxisException
					| PositiveZAxisException e) {
				insideShape += .5;

			}

		}

		double portionOfKnownVolume = insideShape / (double) (coord.size());
		// portion of volume of 3D system
		return portionOfKnownVolume * (2 * a) * (2 * b) * (2 * c);
	}

	/**
	 * Maintenance method for testing
	 * 'generateUniformDistributionOfSamplePoints' method Method prints a file
	 * containing all the generated sample points
	 */
	public void printUniformDistribution() {
		FileWriter file = null;
		try {
			file = new FileWriter("Distribution");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ArrayList<double[]> coord = generateUniformDistributionOfSamplePoints();

		PrintWriter fileOut = new PrintWriter(file);

		for (int i = 0; i < coord.size(); i++) {
			fileOut.println(String.format("[%.4f, %.4f, %.4f]",
					coord.get(i)[x], coord.get(i)[y], coord.get(i)[z]));
		}

	}

	/**
	 * This method can be used on uniformly distributed sample points. For
	 * random distribution use 'determineIfPointIsInsideShapeRandomDistribution'
	 * Checks if point is between specified bounds for theta. Checks if point is
	 * between specified bounds for phi. Checks if point is within the equation
	 * of an ellipsoid. Throws exceptions for cases where the point is on the
	 * edge of the specified shape.
	 * 
	 * fomulas: theta: thetaStart < thetaSample < thetaEnd phi: phiStart <
	 * phiSample < phiEnd equation of ellipsoid: (x^2/a^2) + (y^2/b^2) +
	 * (z^2/c^2) <= 1
	 * 
	 * @param coord
	 * @return
	 */
	private boolean determineIfPointIsInsideShapeNonRandomDistribution(
			double[] coord) throws PointOnEdgeException,
			PositiveZAxisException, NegativeZAxisException {

		boolean insideShape = false;

		try {
			if (pointInsideEquationOfEllipsoidRect(coord)
					&& pointBetweenPhiStartAndPhiEndRect(coord)
					&& pointBetweenThetaStartAndThetaEndEvenDist(coord)) {
				insideShape = true;
			}
		} catch (NegativeZAxisException e) {
			if (radianMeasureOffZAxisEnd == Math.PI) {
				throw new NegativeZAxisException(coord);
			}
		} catch (PositiveZAxisException e1) {
			if (radianMeasureOffZAxisStart == 0) {
				throw new PositiveZAxisException(coord);
			}
		}

		return insideShape;
	}

	/**
	 * checks if theta value of sample point is > thetaStart but < thetaEnd. If
	 * theta is == to startTheta or == to endTheta, a PointOnEdgeException is
	 * thrown
	 * 
	 * @param coord
	 * @return
	 * @throws NegativeZAxisException
	 */
	private boolean pointBetweenThetaStartAndThetaEndEvenDist(double[] coord)
			throws NegativeZAxisException, PositiveZAxisException,
			PointOnEdgeException {
		boolean insideThetaBound = false;

		try {
			if (startRadianTheta <= getThetaValueOfSampleCoordStartTheta(coord)
					&& getThetaValueOfSampleCoordEndTheta(coord) <= endRadianTheta) {

				if (startRadianTheta == getThetaValueOfSampleCoordStartTheta(coord)
						|| getThetaValueOfSampleCoordEndTheta(coord) == endRadianTheta) {
					throw new PointOnEdgeException(coord);
				}

				insideThetaBound = true;
			}
		} catch (PointOnOriginException e) {
			insideThetaBound = true;
			System.out.println(e.getMessage());
		}

		return insideThetaBound;
	}

	/**
	 * Returns an ArrayList<double[]> of 3D sample points of uniform
	 * distribution
	 * 
	 * @return
	 */
	public ArrayList<double[]> generateUniformDistributionOfSamplePoints() {
		final int samplesOnAxis = 50;

		ArrayList<double[]> coord = new ArrayList<double[]>();
		double xVar = 0;
		for (int xAxis = 0; xAxis < samplesOnAxis; xAxis++) {
			double yVar = 0;
			xVar += a / (double) samplesOnAxis;

			for (int yAxis = 0; yAxis < samplesOnAxis; yAxis++) {
				double zVar = 0;
				yVar += b / (double) samplesOnAxis;

				for (int zAxis = 0; zAxis < samplesOnAxis; zAxis++) {

					zVar += c / (double) samplesOnAxis;

					// first octant
					coord.add(new double[] { xVar, yVar, zVar });
					// second octant
					coord.add(new double[] { xVar * -1, yVar, zVar });
					// third octant
					coord.add(new double[] { xVar * -1, yVar * -1, zVar });
					// fourth octant
					coord.add(new double[] { xVar, yVar * -1, zVar });
					// fifth octant
					coord.add(new double[] { xVar, yVar, zVar * -1 });
					// sixth octant
					coord.add(new double[] { xVar * -1, yVar, zVar * -1 });
					// seventh octant
					coord.add(new double[] { xVar * -1, yVar * -1, zVar * -1 });
					// eighth octant
					coord.add(new double[] { xVar, yVar * -1, zVar * -1 });

				}

			}

		}

		// point on origin
		coord.add(new double[] { 0, 0, 0 });
		// add special conditions
		addZAxisPoints(coord, samplesOnAxis);
		addXYPlanePoints(coord, samplesOnAxis);
		addXAxisPoints(coord, samplesOnAxis);
		addYAxisPoints(coord, samplesOnAxis);
		return coord;

	}

	/**
	 * Helper method for uniform distribution method. Adds points above and
	 * below the origin (z axis) to ArrayList
	 * 
	 * @param coord
	 * @param samplesOnAxis
	 */
	private void addZAxisPoints(ArrayList<double[]> coord, int samplesOnAxis) {
		double zVar = 0;
		zVar += c / (double) samplesOnAxis;

		for (int i = 0; i < samplesOnAxis; i++) {
			zVar += c / (double) samplesOnAxis;

			if (i % 2.0 > 0.0) {
				coord.add(new double[] { 0, 0, zVar });
				coord.add(new double[] { 0, 0, -zVar });
			} else {
				coord.add(new double[] { 0, 0, zVar });
				coord.add(new double[] { 0, 0, -zVar });
			}
		}

	}

	/**
	 * Helper method for uniform distribution method. Adds points contained
	 * within the xy plane to the ArrayList
	 * 
	 * @param coord
	 * @param samplesOnAxis
	 */
	private void addXYPlanePoints(ArrayList<double[]> coord, int samplesOnAxis) {
		double xVar = 0;

		for (int i = 0; i < samplesOnAxis; i++) {

			double yVar = 0;

			for (int index = 0; index < samplesOnAxis; index++) {
				yVar += b / (double) samplesOnAxis;

				// first quadrant
				coord.add(new double[] { xVar, yVar, 0 });
				// second quadrant
				coord.add(new double[] { -xVar, yVar, 0 });
				// third quadrant
				coord.add(new double[] { -xVar, -yVar, 0 });
				// fourth quadrant
				coord.add(new double[] { xVar, -yVar, 0 });
			}

		}
	}

	/**
	 * Helper method for uniform distribution method. Adds points above and
	 * below x axis to ArrayList
	 * 
	 * @param coord
	 * @param samplesOnAxis
	 */
	private void addXAxisPoints(ArrayList<double[]> coord, int samplesOnAxis) {
		double xVar = 0;

		for (int i = 0; i < samplesOnAxis; i++) {
			xVar += a / (double) samplesOnAxis;
			double zVar = 0;

			for (int index = 0; index < samplesOnAxis; index++) {
				zVar += c / (double) samplesOnAxis;

				if (index % 2.0 > 0.0) {
					coord.add(new double[] { xVar, 0, zVar });
					coord.add(new double[] { -xVar, 0, -zVar });
				} else {
					coord.add(new double[] { xVar, 0, zVar });
					coord.add(new double[] { -xVar, 0, -zVar });
				}
			}
		}
	}

	/**
	 * Helper method for uniform distribution method. Adds points above and
	 * below x axis to ArrayList
	 * 
	 * @param coord
	 * @param samplesOnAxis
	 */
	private void addYAxisPoints(ArrayList<double[]> coord, int samplesOnAxis) {
		double yVar = 0;

		for (int i = 0; i < samplesOnAxis; i++) {
			yVar += b / (double) samplesOnAxis;
			double zVar = 0;

			for (int index = 0; index < samplesOnAxis; index++) {
				zVar += c / (double) samplesOnAxis;

				coord.add(new double[] { 0, yVar, zVar });
				coord.add(new double[] { 0, -yVar, -zVar });
			}
		}
	}

	// ---------------------------Tiny Volumes---------------------------
	/**
	 * Uses a uniformly distributed set of 3 dimensional coordinates to estimate
	 * volume. Formula for number of points is = [(number of sample points on 1
	 * positive axis)^3] * 8. Due to computation time, samples are limited to
	 * ~2,700,000 sample points.
	 * 
	 * @return
	 */
	public double getEstimatedVolumeThruUniformDistributionTinyVolumesSum() {

		ArrayList<double[]> coord = generateUniformDistributionOfSamplePoints();

		double rectangularPrisimVolume = this.a / this.sampleSize * this.b
				/ this.sampleSize * this.c / this.sampleSize;
		double shapeVolume = 0;

		for (int i = 0; i < coord.size(); i++) {
			try {
				if (determineIfPointIsInsideShapeNonRandomDistribution(coord
						.get(i))) {
					shapeVolume += rectangularPrisimVolume;
				}
			} catch (PointOnEdgeException | PositiveZAxisException e) {
				shapeVolume += rectangularPrisimVolume / 2;
			}

		}

		// portion of volume of 3D system
		return shapeVolume;
	}

}
