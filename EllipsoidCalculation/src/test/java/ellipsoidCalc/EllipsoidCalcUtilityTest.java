package ellipsoidCalc;

import static org.junit.Assert.*;
import static org.hamcrest.CoreMatchers.*;

import org.junit.Test;

public class EllipsoidCalcUtilityTest {

	@Test
	public void testConversionToDecimal() {
		EllipsoidCalcUtility ellipse = new EllipsoidCalcUtility();
		assertThat(ellipse.convertToDecimal("3/4"), is(.75));
	}
	
	@Test
	public void testConversionToRadians(){
		EllipsoidCalcUtility ellipse = new EllipsoidCalcUtility();
		assertThat(ellipse.convertThetaToRadians("180", true), is(Math.PI));
	}
	
	@Test
	public void testVolumeEighth(){
		EllipsoidCalcUtility ellipse = new EllipsoidCalcUtility();
		double quarterOfEllipseVolume = (1/(double)6) * Math.PI * 2 * 2 * 3;
		assertThat(ellipse.getVolume(0, (Math.PI * .5), (Math.PI * .5), 0, 2, 2, 3), is(quarterOfEllipseVolume));
	}
	
	@Test
	public void testVolumeQuarter(){
		EllipsoidCalcUtility ellipse = new EllipsoidCalcUtility();
		double quarterOfEllipseVolume = (1/(double)3) * Math.PI * 2 * 2 * 3;
		assertThat(ellipse.getVolume(0, (Math.PI * .5), Math.PI, 0, 2, 2, 3), is(quarterOfEllipseVolume));
	}
	
	@Test
	public void testVolumeHalf(){
		EllipsoidCalcUtility ellipse = new EllipsoidCalcUtility();
		double halfOfEllipseVolume = (2/(double)3) * Math.PI * 2 * 2 * 3;
		assertThat(ellipse.getVolume(0, (Math.PI), Math.PI, 0, 2, 2, 3), is(halfOfEllipseVolume));
	}
	
	@Test
	public void testVolumeWhole(){
		EllipsoidCalcUtility ellipse = new EllipsoidCalcUtility();
		double halfOfEllipseVolume = (4/(double)3) * Math.PI * 2 * 2 * 3;
		assertThat(ellipse.getVolume(0, (Math.PI * 2), Math.PI, 0, 2, 2, 3), is(halfOfEllipseVolume));
	}

}
