package ellipsoidCalc;

import static org.junit.Assert.*;
import static org.hamcrest.CoreMatchers.*;

import org.junit.Test;

public class EllipsoidCalcUtilityTest {

	@Test
	public void test() {
		EllipsoidCalcUtility ellipse = new EllipsoidCalcUtility();
		assertThat(ellipse.convertToDecimal("3/4"), is(.75));
	}

}
