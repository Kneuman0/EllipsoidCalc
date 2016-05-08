package biz.personalAcademics.ellipsoidCalc;

@SuppressWarnings("serial")
public class PositiveZAxisException extends Exception{
	
	public PositiveZAxisException(double[] coord) {
		super(String.format("Point on negative z Axis\nx = %.4f, y = %.4f, z = %.4f", coord[0], coord[1], coord[2]));
	}

}
