package biz.personalAcademics.ellipsoid.customExceptions;

@SuppressWarnings("serial")
public class PointOnOriginException extends RuntimeException {
	public PointOnOriginException(double[] coord) {
		super(String.format("Origin error for point\nx = %.4f, y = %.4f, z = %.4f", coord[0], coord[1], coord[2]));
	}
}
