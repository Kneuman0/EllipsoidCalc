package biz.personalAcademics.ellipsoidCalc;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;

public class Ellipsoid {

	private double startRadianTheta, endRadianTheta, radianMeasureOffZAxisEnd,
			radianMeasureOffZAxisStart, a, b, c;

	private int x, y, z;

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
		this.x = 0;
		this.y = 1;
		this.z = 2;

		/*
		 * the nextDouble method will be called from this object and multiplied
		 * by one of the axes. This will yield a coordinate that lays within the
		 * 1st octant.
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

	/**
	 * Estimates the volume of the ellipsoid using the specified number of
	 * sample points
	 * 
	 * @param sampleSize
	 * @return
	 */
	public double getEstimatedVolume(int sampleSize) {
		int insideShape = 0;

		for (int i = 0; i < sampleSize; i++) {
			if (determineIfPointIsInsideShape(generateRandom3DSamplePoint())) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) sampleSize;

		// portion of volume of 3D system
		return portionOfKnownVolume * (2 * a) * (2 * b) * (2 * c);
	}

	/**
	 * Estimates the volume of the portion of the ellipsoid using 500000 sample
	 * points
	 * 
	 * @return
	 */
	public double getEstimatedVolume() {
		int insideShape = 0;
		final int sampleSize = 5_000_000;

		for (int i = 0; i < sampleSize; i++) {
			if (determineIfPointIsInsideShape(generateRandom3DSamplePoint())) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) sampleSize;

		// portion of volume of 3D system
		return portionOfKnownVolume * (2 * a) * (2 * b) * (2 * c);
	}
	
	public double getEstimatedVolumeThruUniformDistribution(){
		int insideShape = 0;
		ArrayList<double[]> coord = generateUniformDistributionOfSamplePoints();

		for (int i = 0; i < coord.size(); i++) {
			if (determineIfPointIsInsideShape(coord.get(i))) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) (coord.size());
		System.out.println(portionOfKnownVolume);
		System.out.println(insideShape);

		// portion of volume of 3D system
		return portionOfKnownVolume * (2 * a) * (2 * b) * (2 * c);
	}

	private boolean determineIfPointIsInsideShape(double[] coord) {
		boolean insideShape = false;

		try {
			if (pointInsideEquationOfEllipsoid(coord)
					&& pointBetweenPhiStartAndPhiEnd(coord)
					&& pointBetweenThetaStartAndThetaEnd(coord)) {
				insideShape = true;
			}
		} catch (NegativeZAxisException e) {
			if (radianMeasureOffZAxisEnd == Math.PI) {
				insideShape = true;
			}
		} catch (PositiveZAxisException e1) {
			if (radianMeasureOffZAxisStart == Math.PI / 2) {
				insideShape = true;
			}
		}

		return insideShape;
	}

	private boolean pointInsideEquationOfEllipsoid(double[] coord) {
		boolean insideShape = false;

		// point is inside shape if x^2/a^2 + y^2/b^2 + z^2/c^2 <= 1
		double valueOfPointInEllipsoidEquation = (coord[x] * coord[x])
				/ (a * a) + (coord[y] * coord[y]) / (b * b)
				+ (coord[z] * coord[z]) / (c * c);

		if (valueOfPointInEllipsoidEquation <= 1) {
			insideShape = true;
		}

		return insideShape;
	}

	private boolean pointBetweenPhiStartAndPhiEnd(double[] coord) {
		boolean insidePhiBound = false;

		try {
			if (radianMeasureOffZAxisStart <= getPhiValueOfSampleCoord(coord)
					&& getPhiValueOfSampleCoord(coord) <= radianMeasureOffZAxisEnd) {
				insidePhiBound = true;
			}
		} catch (PointOnOriginException e) {
			insidePhiBound = true;
			System.out.println(e.getMessage());
		}
		return insidePhiBound;
	}

	private boolean pointBetweenThetaStartAndThetaEnd(double[] coord)
			throws NegativeZAxisException {
		boolean insideThetaBound = false;

		try {
			if (startRadianTheta <= getThetaValueOfSampleCoordStartTheta(coord)
					&& getThetaValueOfSampleCoordEndTheta(coord) <= endRadianTheta) {
				insideThetaBound = true;
			}
		} catch (PointOnOriginException e) {
			insideThetaBound = true;
			System.out.println(e.getMessage());
		}

		return insideThetaBound;
	}

	private double getPhiValueOfSampleCoord(double[] coord) {
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
					+ Math.PI;

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
					+ Math.PI;

			// in zx plane above xy plane
		} else if (coord[z] > 0 && coord[x] != 0 && coord[y] == 0) {

			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])));

			// in zx plane below xy plane
		} else if (coord[z] < 0 && coord[x] != 0 && coord[y] == 0) {

			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])))
					+ Math.PI;

			// in xy plane
		} else if (coord[z] == 0) {
			return Math.PI / 2;

			// on origin
		} else {

			throw new PointOnOriginException(coord);
		}

		// do case where phi is below x y plane

	}

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

			// on origin
		} else if (coord[x] == 0 && coord[y] == 0 && coord[z] < 0) {

			throw new NegativeZAxisException(coord);

		} else if (coord[x] == 0 && coord[y] == 0 && coord[z] > 0) {

			throw new PositiveZAxisException(coord);

			// on origin
		} else {
			throw new PointOnOriginException(coord);
		}
	}

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

		} else if (coord[x] < 0 && coord[y] == 0 && coord[z] < 0) {

			throw new NegativeZAxisException(coord);

		} else if (coord[x] < 0 && coord[y] == 0 && coord[z] > 0) {

			throw new PositiveZAxisException(coord);

			// on origin
		} else {
			throw new PointOnOriginException(coord);
		}
	}

	private int getRandomNegation() {
		boolean negativeSign = portionOfAxis.nextBoolean();

		if (negativeSign) {
			return 1;
		} else {
			return -1;
		}
	}

	/**
	 * generates a 3 dimensional sample point that lays with the first octant of
	 * the ellipsoid Cartesian field.
	 * 
	 * @return
	 */
	private double[] generateRandom3DSamplePoint() {
		double[] samplePoint = new double[3];

		samplePoint[x] = a * portionOfAxis.nextDouble() * getRandomNegation();
		samplePoint[y] = b * portionOfAxis.nextDouble() * getRandomNegation();
		samplePoint[z] = c * portionOfAxis.nextDouble() * getRandomNegation();

		return samplePoint;
	}

	public ArrayList<double[]> generateUniformDistributionOfSamplePoints() {
		final int samplesOnAxis = 150;

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

	private void addZAxisPoints(ArrayList<double[]> coord, int samplesOnAxis) {
		double zVar = 0;
		zVar += c / (double) samplesOnAxis;

		for (int i = 0; i < samplesOnAxis; i++) {
			zVar += c / (double) samplesOnAxis;
			coord.add(new double[] { 0, 0, zVar });
			coord.add(new double[] { 0, 0, -zVar });
		}

	}

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

	private void addXAxisPoints(ArrayList<double[]> coord, int samplesOnAxis) {
		double xVar = 0;

		for (int i = 0; i < samplesOnAxis; i++) {
			xVar += a / (double) samplesOnAxis;
			double zVar = 0;

			for (int index = 0; index < samplesOnAxis; index++) {
				zVar += c / (double) samplesOnAxis;
				coord.add(new double[] { xVar, 0, zVar });
				coord.add(new double[] { -xVar, 0, -zVar });
			}
		}
	}
	
	private void addYAxisPoints(ArrayList<double[]> coord, int samplesOnAxis){
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
	
	public void printUniformDistribution(){
		FileWriter file = null;
		try {
			file = new FileWriter("Distribution");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ArrayList<double[]> coord = generateUniformDistributionOfSamplePoints();
		
		PrintWriter fileOut = new PrintWriter(file);
		
		for(int i = 0; i < coord.size(); i++){
			fileOut.println(String.format("[%.4f, %.4f, %.4f]", 
					coord.get(i)[x], coord.get(i)[y], coord.get(i)[z]));
		}
		
	}

}
