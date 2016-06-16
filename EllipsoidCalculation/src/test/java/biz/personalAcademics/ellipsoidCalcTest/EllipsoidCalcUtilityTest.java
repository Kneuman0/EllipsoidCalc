package biz.personalAcademics.ellipsoidCalcTest;

import static org.junit.Assert.*;
import static org.hamcrest.CoreMatchers.*;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import biz.personalAcademics.ellipsoid.Ellipsoid;
import biz.personalAcademics.ellipsoid.EllipsoidCylindricalCoords;
import biz.personalAcademics.ellipsoid.EllipsoidEvenDistributionRectCoords;
import biz.personalAcademics.ellipsoid.EllipsoidRectangularCoords;
import biz.personalAcademics.ellipsoid.EllipsoidSphericalCoords;
import biz.personalAcademics.ellipsoid.customExceptions.InvalidUserInputException;

public class EllipsoidCalcUtilityTest {
	
	private final boolean testRectangularRandomDist = false;
	private final boolean testSphericalRandomDist = true;
	private final boolean testMonteCarloUniformDist = false;
	private final boolean testCylindricalRandomDist = false;
	private final double A =  100, B = 8, C = 8;

	@Rule
	public ExpectedException invalidInput = ExpectedException.none();

//---------------------Random Dis. Cylinder-----------------------------
	@Test
	public void testEstimatedEighthVolumeCalculationRandomCylinder(){
		
		if(testCylindricalRandomDist){
		
			double EighthOfEllipseVolume = (1 / (double) 6) * Math.PI * A * B * C;
			EllipsoidCylindricalCoords ellip = new EllipsoidCylindricalCoords(0, (.5 * Math.PI), 0, (Math.PI/2),
					A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", EighthOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}
	
	/**
	 * Using formula for ellipsoid (1/4) * (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedQuarterVolumeCalculationRandomCylinder() {
		
		if(testCylindricalRandomDist){
			double quarterOfEllipseVolume = (1 / (double) 3) * Math.PI * A * B * C;
			EllipsoidCylindricalCoords ellip = new EllipsoidCylindricalCoords((0), (Math.PI/2), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", quarterOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}

	/**
	 * Using formula for ellipsoid (1/2) * (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedHalfVolumeCalculationRandomCylinder() {
		
		if(testCylindricalRandomDist){
			double halfOfEllipseVolume = (2 / (double) 3) * Math.PI * A * B * C;
			EllipsoidCylindricalCoords ellip = new EllipsoidCylindricalCoords(0, (Math.PI), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", halfOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}

	/**
	 * Using formula for ellipsoid (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedWholehVolumeCalculationRandomCylinder() {
		
		if(testCylindricalRandomDist){
			double fullEllipseVolume = (4 / (double) 3) * Math.PI * A * B * C;
			EllipsoidCylindricalCoords ellip = new EllipsoidCylindricalCoords(0, (Math.PI * 2), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", fullEllipseVolume)));
			
		}else{
			// if not testing, test returns true
			assert(true);
			
		}	
	}
	
//---------------------Random Dist. Sphere--------------------------------
	
	@Test
	public void testEstimatedEighthVolumeCalculationRandomSphere(){
		
		if(testSphericalRandomDist){
		
			double EighthOfEllipseVolume = (1 / (double) 6) * Math.PI * A * B * C;
			EllipsoidSphericalCoords ellip = new EllipsoidSphericalCoords(0, (.5 * Math.PI), 0, (Math.PI/2),
					A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", EighthOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}
	
	/**
	 * Using formula for ellipsoid (1/4) * (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedQuarterVolumeCalculationRandomSphere() {
		
		if(testSphericalRandomDist){
			double quarterOfEllipseVolume = (1 / (double) 3) * Math.PI * A * B * C;
			EllipsoidSphericalCoords ellip = new EllipsoidSphericalCoords((0), (Math.PI/2), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", quarterOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}

	/**
	 * Using formula for ellipsoid (1/2) * (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedHalfVolumeCalculationRandomSphere() {
		
		if(testSphericalRandomDist){
			double halfOfEllipseVolume = (2 / (double) 3) * Math.PI * A * B * C;
			EllipsoidSphericalCoords ellip = new EllipsoidSphericalCoords(0, (Math.PI), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", halfOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}

	/**
	 * Using formula for ellipsoid (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedWholehVolumeCalculationRandomSphere() {
		
		if(testSphericalRandomDist){
			double fullEllipseVolume = (4 / (double) 3) * Math.PI * A * B * C;
			EllipsoidSphericalCoords ellip = new EllipsoidSphericalCoords(0, (Math.PI * 2), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", fullEllipseVolume)));
			
		}else{
			// if not testing, test returns true
			assert(true);
			
		}	
	}
	
//---------------------Random Dist. Rectangular Prism-----------------------
	@Test
	public void testEstimatedEighthVolumeCalculationRandomRect(){
		
		if(testRectangularRandomDist){
			double EighthOfEllipseVolume = (1 / (double) 6) * Math.PI * A * B * C;
			EllipsoidRectangularCoords ellip = new EllipsoidRectangularCoords(0, (.5 * Math.PI), 0, (Math.PI/2),
					A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", EighthOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}	
	}
	
	/**
	 * Using formula for ellipsoid (1/4) * (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedQuarterVolumeCalculationRandomRect() {
		
		if(testRectangularRandomDist){
			double quarterOfEllipseVolume = (1 / (double) 3) * Math.PI * A * B * C;
			EllipsoidRectangularCoords ellip = new EllipsoidRectangularCoords((0), (Math.PI/2), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", quarterOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}	
	}

	/**
	 * Using formula for ellipsoid (1/2) * (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedHalfVolumeCalculationRandomRect() {
		
		if(testRectangularRandomDist){
			double halfOfEllipseVolume = (2 / (double) 3) * Math.PI * A * B * C;
			EllipsoidRectangularCoords ellip = new EllipsoidRectangularCoords(0, (Math.PI), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", halfOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}	
	}

	/**
	 * Using formula for ellipsoid (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedWholehVolumeCalculationRandomRect() {
		
		if(testRectangularRandomDist){
			double fullEllipseVolume = (4 / (double) 3) * Math.PI * A * B * C;
			EllipsoidRectangularCoords ellip = new EllipsoidRectangularCoords(0, (Math.PI * 2), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", fullEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}
	
//------------------Uniform (equally spaced) dist------------------
	@Test
	public void testEstimatedEighthVolumeCalculationUniformDist(){
		
		if(testMonteCarloUniformDist){
			double EighthOfEllipseVolume = (1 / (double) 6) * Math.PI * A * B * C;
			EllipsoidEvenDistributionRectCoords ellip = new EllipsoidEvenDistributionRectCoords(0, (.5 * Math.PI), 0, (Math.PI/2),
					A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", EighthOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}
	
	/**
	 * Using formula for ellipsoid (1/4) * (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedQuarterVolumeCalculationUniformDist() {
		
		if(testMonteCarloUniformDist){
			double quarterOfEllipseVolume = (1 / (double) 3) * Math.PI * A * B * C;
			EllipsoidEvenDistributionRectCoords ellip = new EllipsoidEvenDistributionRectCoords((0), (Math.PI/2), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", quarterOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}

	/**
	 * Using formula for ellipsoid (1/2) * (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedHalfVolumeCalculationUniformDist() {
		
		if(testMonteCarloUniformDist){
			double halfOfEllipseVolume = (2 / (double) 3) * Math.PI * A * B * C;
			EllipsoidEvenDistributionRectCoords ellip = new EllipsoidEvenDistributionRectCoords(0, (Math.PI), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", halfOfEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}

	/**
	 * Using formula for ellipsoid (4/3)pi(abc)
	 */
	@Test
	public void testEstimatedWholehVolumeCalculationUniformDist() {
		
		if(testMonteCarloUniformDist){
			double fullEllipseVolume = (4 / (double) 3) * Math.PI * A * B * C;
			EllipsoidEvenDistributionRectCoords ellip = new EllipsoidEvenDistributionRectCoords(0, (Math.PI * 2), 0, Math.PI, A, B, C);
		
			assertThat(String.format("%.2f", ellip.getEstimatedVolume()),
					is(String.format("%.2f", fullEllipseVolume)));
		}else{
			// if not testing, test returns true
			assert(true);
			
		}
	}
	
//------------------------Exact Volume---------------------------
	
	/**
	 * Using formula for ellipsoid (1/8) * (4/3)pi(abc)
	 */
	@Test
	public void testVolumeEighth() {
		double EighthOfEllipseVolume = (1 / (double) 6) * Math.PI * A * B * C;
		Ellipsoid ellip = new Ellipsoid(0, (.5 * Math.PI), 0, (Math.PI/2),
				A, B, C);
		
		assertThat(String.format("%.13f", ellip.getEstimatedVolume()),
				is(String.format("%.13f", EighthOfEllipseVolume)));
	}

	/**
	 * Using formula for ellipsoid (1/4) * (4/3)pi(abc)
	 */
	@Test
	public void testVolumeQuarter() {
		double quarterOfEllipseVolume = (1 / (double) 3) * Math.PI * A * B * C;
		Ellipsoid ellip = new Ellipsoid((0), (Math.PI/2), 0, Math.PI, A, B, C);
		
		assertThat(String.format("%.13f", ellip.getEstimatedVolume()),
				is(String.format("%.13f", quarterOfEllipseVolume)));
	}

	/**
	 * Using formula for ellipsoid (1/2) * (4/3)pi(abc)
	 */
	@Test
	public void testVolumeHalf() {
		double halfOfEllipseVolume = (2 / (double) 3) * Math.PI * A * B * C;
		Ellipsoid ellip = new Ellipsoid(0, (Math.PI), 0, Math.PI, A, B, C);
		
		assertThat(String.format("%.13f", ellip.getEstimatedVolume()),
				is(String.format("%.13f", halfOfEllipseVolume)));
	}

	/**
	 * Using formula for ellipsoid (4/3)pi(abc)
	 */
	@Test
	public void testVolumeWhole() {
		double fullEllipseVolume = (4 / (double) 3) * Math.PI * A * B * C;
		Ellipsoid ellip = new Ellipsoid(0, (Math.PI * 2), 0, Math.PI, A, B, C);
		
		assertThat(String.format("%.13f", ellip.getEstimatedVolume()),
				is(String.format("%.13f", fullEllipseVolume)));
	}
	
//--------------------------Other tests---------------------------
	
	
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

}
