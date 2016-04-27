package biz.personalAcademics.ellipsoidCalc;

public class NegativeZAxisException extends RuntimeException{
		
	public NegativeZAxisException(double[] coord) {
		super(String.format("Point on negative z Axis\nx = %.4f, y = %.4f, z = %.4f", coord[0], coord[1], coord[2]));
	}
}
