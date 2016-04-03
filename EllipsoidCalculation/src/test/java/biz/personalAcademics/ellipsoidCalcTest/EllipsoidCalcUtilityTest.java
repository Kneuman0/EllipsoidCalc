package biz.personalAcademics.ellipsoidCalcTest;

import static org.junit.Assert.*;
import static org.hamcrest.CoreMatchers.*;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import biz.personalAcademics.ellipsoidCalc.*;

public class EllipsoidCalcUtilityTest {

	@Rule
	public ExpectedException invalidInput = ExpectedException.none();

	@Test
	public void testConversionToDecimal() {
		assertThat(Ellipsoid.convertToDecimal("3/4"), is(.75));
	}

	@Test
	public void testConversionToRadiansCalculation() {
		assertThat(Ellipsoid.convertThetaToRadians("180", true), is(Math.PI));
	}

	@Test
	public void testConversionToRadiansForInvalidUserInput() {
		invalidInput.expect(InvalidUserInputException.class);
		invalidInput.expectMessage(containsString("abcd"));
		Ellipsoid.convertThetaToRadians("abcd", true);
	}

	/**
	 * Using formula for ellipsoid (1/8) * (4/3)pi(abc)
	 */
	@Test
	public void testVolumeEighth() {
		double EighthOfEllipseVolume = (1 / (double) 6) * Math.PI * 2 * 2 * 3;
		Ellipsoid ellip = new Ellipsoid(0, (Math.PI * .5), (Math.PI * .5), 0,
				2, 2, 3);
		
		assertThat(String.format("%.13f", ellip.getVolume()),
				is(String.format("%.13f", EighthOfEllipseVolume)));
	}

	/**
	 * Using formula for ellipsoid (1/4) * (4/3)pi(abc)
	 */
	@Test
	public void testVolumeQuarter() {
		double quarterOfEllipseVolume = (1 / (double) 3) * Math.PI * 2 * 2 * 3;
		Ellipsoid ellip = new Ellipsoid(0, (Math.PI * .5), Math.PI, 0, 2, 2, 3);
		
		assertThat(String.format("%.13f", ellip.getVolume()),
				is(String.format("%.13f", quarterOfEllipseVolume)));
	}

	/**
	 * Using formula for ellipsoid (1/2) * (4/3)pi(abc)
	 */
	@Test
	public void testVolumeHalf() {
		double halfOfEllipseVolume = (2 / (double) 3) * Math.PI * 2 * 2 * 3;
		Ellipsoid ellip = new Ellipsoid(0, (Math.PI), Math.PI, 0, 2, 2, 3);
		
		assertThat(String.format("%.13f", ellip.getVolume()),
				is(String.format("%.13f", halfOfEllipseVolume)));
	}

	/**
	 * Using formula for ellipsoid (4/3)pi(abc)
	 */
	@Test
	public void testVolumeWhole() {
		double fullEllipseVolume = (4 / (double) 3) * Math.PI * 2 * 2 * 3;
		Ellipsoid ellip = new Ellipsoid(0, (Math.PI * 2), Math.PI, 0, 2, 2, 3);
		
		assertThat(String.format("%.13f", ellip.getVolume()),
				is(String.format("%.13f", fullEllipseVolume)));
	}

}
