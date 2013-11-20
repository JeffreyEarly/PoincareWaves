//
//  main.m
//  PoincareWaves
//
//  Created by Jeffrey J. Early on 11/19/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

int main(int argc, const char * argv[])
{

	@autoreleasepool {
		// Reasonable parameters to nondimensionalize by.
		GLFloat U_max = 0.01;                // m/s
		GLFloat L_domain = 4000;            // m
		GLFloat L_r = 100;                  // m
		GLFloat latitude = 33;              // degrees
        GLFloat maxInertialPeriods = 5;     // # of inertial periods
        GLFloat sampleTimeInMinutes = 30;  // output sample time.
		
		GLFloat f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
		GLFloat g = 9.81;
		GLFloat D = pow(L_r*f0,2)/g;
		GLFloat c = sqrt(g*D);
        
        // This is good for unit testing.
        BOOL shouldEvolveInertialOnly = NO;
		
		/************************************************************************************************/
		/*		Define the problem dimensions															*/
		/************************************************************************************************/
		
		GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:512 domainMin: -L_domain/2 length:L_domain];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:512 domainMin: -L_domain/2  length:L_domain];
		yDim.name = @"y";
		GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
		
		// Variables are always tied to a particular equation---so we create an equation object first.
		GLEquation *equation = [[GLEquation alloc] init];
		
		NSArray *spatialDimensions = @[xDim, yDim];
		GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
		GLFunction *y = [GLFunction functionOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
		
		/************************************************************************************************/
		/*		Create and cache the differential operators we will need								*/
		/************************************************************************************************/
		
		NSArray *spectralDimensions = [x dimensionsTransformedToBasis: x.differentiationBasis];
		GLDimension *kDim = spectralDimensions[0];
		GLDimension *lDim = spectralDimensions[1];
		GLFunction *k = [GLFunction functionOfRealTypeFromDimension: kDim withDimensions: spectralDimensions forEquation: equation];
		GLFunction *l = [GLFunction functionOfRealTypeFromDimension: lDim withDimensions: spectralDimensions forEquation: equation];
		
		GLFunction *kH = [[[[k multiply: k] plus: [l multiply: l]] sqrt] times: @(2*M_PI)];
	    GLFunction *omega = [[[[kH multiply: kH] multiply: @(g*D)] plus: @(f0*f0)] sqrt];
		
		GLFloat minPeriod = 2*M_PI/[omega maxNow];
		GLFloat maxPeriod = 2*M_PI/[omega minNow];
		
		NSLog(@"The minimum wave period is %.0f minutes (based on the horizontal resolution), maximum is %.1f hours (based on the Coriolis frequency).", round(minPeriod/60), maxPeriod/(60*60));
		NSLog(@"Sampling for longer than %.1f days will expose the time discretization (this is the domainWidth/c).", L_domain/(c*86400));
        NSLog(@"You have chosen to sample at a rate of %.0f minutes, for %.1f days", sampleTimeInMinutes, maxInertialPeriods*2*M_PI/(f0*86400));
        NSLog(@"Your inertial circles will have a radius of %.0f meters", U_max/f0);
		
		/************************************************************************************************/
		/*		Now set the magnitude of the wave at each component.									*/
		/************************************************************************************************/
		
		// Compute the angle of the wavevector...
		GLFunction *alpha = [l atan2: k];
		
		// Set the magnitude of the components
		GLFunction *U_mag = [[[omega multiply: @(1/f0)] pow: -2.5] multiply: @(2*U_max)];
		if (shouldEvolveInertialOnly) {
            [U_mag zero];
            U_mag = [U_mag setValue: 2*U_max atIndices: @"0,0"];
        }
        
//        GLFunction *filter = [GLBasisTransformOperation antialiasingFilterFromDimensions: spectralDimensions forEquation: equation];
//        U_mag = [U_mag multiply: filter];
		
		// This computes how they map to u,v
		GLFunction *u_phase = [[alpha cos] minus: [[[[alpha sin] dividedBy: omega] multiply: @(f0)] swapComplex]];
		GLFunction *v_phase = [[alpha sin] plus: [[[[alpha cos] dividedBy: omega] multiply: @(f0)] swapComplex]];
		
		u_phase = [U_mag multiply: u_phase];
		v_phase = [U_mag multiply: v_phase];
		
		// Initial randomized phase
		GLFunction *phi0 = [GLFunction functionWithRandomValuesBetween: -M_PI and: M_PI withDimensions: spectralDimensions forEquation: equation];
		if (shouldEvolveInertialOnly) [phi0 zero];
		
		/************************************************************************************************/
		/*		Now run the model																		*/
		/************************************************************************************************/
		
		GLScalar *t = [GLScalar scalarWithValue: 0.75*2*M_PI/f0 forEquation: equation];
		GLFunction *phi = [[omega multiply: t] plus: phi0];
		GLFunction *time_phase = [[[phi cos] multiply: @(0.5)] plus: [[[phi sin] multiply: @(0.5)] swapComplex]];
		
		GLFunction *u = [[u_phase multiply: time_phase] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis)]];
		GLFunction *v = [[v_phase multiply: time_phase] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis)]];
        
        /************************************************************************************************/
		/*		Let's also plop a float at each grid point.                                             */
		/************************************************************************************************/
        
        GLDimension *xFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:xDim.nPoints/8 domainMin: -L_domain/2 length:L_domain];
		xFloatDim.name = @"x-float";
		GLDimension *yFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:yDim.nPoints/8 domainMin: -L_domain/2  length:L_domain];
		yFloatDim.name = @"y-float";
        
        NSArray *floatDimensions = @[xFloatDim, yFloatDim];
		GLFunction *xFloat = [GLFunction functionOfRealTypeFromDimension: xFloatDim withDimensions: floatDimensions forEquation: equation];
		GLFunction *yFloat = [GLFunction functionOfRealTypeFromDimension: yFloatDim withDimensions: floatDimensions forEquation: equation];
        
		GLFunction *xPosition = [GLFunction functionFromFunction: xFloat];
		GLFunction *yPosition = [GLFunction functionFromFunction: yFloat];
        
		/************************************************************************************************/
		/*		Create a NetCDF file and mutable variables in order to record some of the time steps.	*/
		/************************************************************************************************/
		
		NSString *path = [[NSSearchPathForDirectoriesInDomains(NSDesktopDirectory, NSUserDomainMask, YES) objectAtIndex:0] stringByAppendingPathComponent:@"PoincareWaves.nc"];
		GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: [NSURL URLWithString: path] forEquation: equation overwriteExisting: YES];
		
        [netcdfFile setGlobalAttribute: @(U_max) forKey: @"U_max"];
		[netcdfFile setGlobalAttribute: @(L_domain) forKey: @"L_domain"];
		[netcdfFile setGlobalAttribute: @(L_r) forKey: @"L_r"];
        [netcdfFile setGlobalAttribute: @(latitude) forKey: @"latitude"];
        [netcdfFile setGlobalAttribute: @(f0) forKey: @"f0"];
        [netcdfFile setGlobalAttribute: @(D) forKey: @"D"];
        
		GLMutableVariable *uHistory = [u variableByAddingDimension: tDim];
		uHistory.name = @"u";
		uHistory = [netcdfFile addVariable: uHistory];
		
		GLMutableVariable *vHistory = [v variableByAddingDimension: tDim];
		vHistory.name = @"v";
		vHistory = [netcdfFile addVariable: vHistory];
        
        GLMutableVariable *xPositionHistory = [xPosition variableByAddingDimension: tDim];
		xPositionHistory.name = @"x-position";
		xPositionHistory = [netcdfFile addVariable: xPositionHistory];
        
		GLMutableVariable *yPositionHistory = [yPosition variableByAddingDimension: tDim];
		yPositionHistory.name = @"y-position";
		yPositionHistory = [netcdfFile addVariable: yPositionHistory];
        
        /************************************************************************************************/
		/*		Create the integration object.															*/
        /*      The wave equation is solved analytically, but the particles need to be integrated.      */
		/************************************************************************************************/
        
        CGFloat cfl = 0.5;
        GLFloat timeStep = cfl * xDim.sampleInterval / U_max;
        GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: @[xPosition, yPosition] stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
            
            GLFunction *phi = [[omega multiply: time] plus: phi0];
            GLFunction *time_phase = [[[phi cos] multiply: @(0.5)] plus: [[[phi sin] multiply: @(0.5)] swapComplex]];
            
            GLFunction *u = [[u_phase multiply: time_phase] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis)]];
            GLFunction *v = [[v_phase multiply: time_phase] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis)]];
        
			GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[u, v] secondOperand: yNew];
			
			return interp.result;
		}];
        
        GLFloat maxTime = maxInertialPeriods*2*M_PI/f0;
        GLFloat sampleTime = sampleTimeInMinutes*60;
        for (GLFloat time = sampleTime; time < maxTime; time += sampleTime)
        {
            @autoreleasepool {
                NSLog(@"Logging day: %f, step size: %f.", (time/86400), integrator.lastStepSize/86400);
                
                NSArray *yout = [integrator stepForwardToTime: time];
                
                t = [GLScalar scalarWithValue: time forEquation: equation];
                phi = [[omega multiply: t] plus: phi0];
                time_phase = [[[phi cos] multiply: @(0.5)] plus: [[[phi sin] multiply: @(0.5)] swapComplex]];
                
                u = [[u_phase multiply: time_phase] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis)]];
                v = [[v_phase multiply: time_phase] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis)]];
                
				[tDim addPoint: @(time)];
				[uHistory concatenateWithLowerDimensionalVariable: u alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [vHistory concatenateWithLowerDimensionalVariable: v alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [xPositionHistory concatenateWithLowerDimensionalVariable: yout[0] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                [yPositionHistory concatenateWithLowerDimensionalVariable: yout[1] alongDimensionAtIndex:0 toIndex: (tDim.nPoints-1)];
                
                [netcdfFile waitUntilAllOperationsAreFinished];
            }
        }
		
		[netcdfFile close];
	}
    return 0;
}

