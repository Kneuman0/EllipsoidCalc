package biz.personalAcademics.ellipsoidCalc;

@SuppressWarnings("serial")
public class PointOnEdgeException extends RuntimeException{
	
	public PointOnEdgeException(double[] coord){
		super(String.format("Point on edge\nx = %.f, y = %.f, z = %.f", coord[0], coord[1], coord[2]));
	}
}
