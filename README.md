# EllipsoidCalc

README

Fully functioning stand-alone program (once exported to a runnable .jar file). This program calculates the volume of an ellipsoid (egg-shape) using the length of the axes and spherical coordinates as inputs. Knowledge of spherical coordinates is required to use the program. The Application uses a Monte Carlo algorithm to estimate the volume of the shape in most cases. The sample size is 10 million points. 

The algorithm is limited to 99.0% accuracy when the smallest axis is >=1 and 99.95% accuracy when the smallest axis is >= 10. This is limited due to the compromise between speed of the program and accuracy of the results. 

If the bounds of integration for φ and θ are increments of π/2, the calculation is very quick and accurate. If the axes are as follows:  a = b = c (a sphere is entered), the calculation is also very quick and accurate.

To run the application, navigate to the class named EllipsoidCalcMain and run it.

If you'd like to review the math, navigate to the Documentation and Notes on Math folder and click on: Summary of EllipsoidCalc




# Further reading

	1. Hand Written Formula and Process Volume Calc (legacy document)

