package biz.personalAcademics.debug;

import java.util.Arrays;

import biz.personalAcademics.ellipsoid.Ellipsoid;
import static java.lang.Math.PI;;

public class TestMain {

	public static void main(String[] args) {
		Ellipsoid ellip = new Ellipsoid(0, .5 * PI, 0, .5 * PI, 7, 5, 3);
		System.out.println(Arrays.toString(ellip.getSortedAxes()));
		ellip.setB(3);
		System.out.println(Arrays.toString(ellip.getSortedAxes()));
	}

}
