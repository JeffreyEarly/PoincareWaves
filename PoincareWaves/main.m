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
        BOOL shouldUnitTest = YES;
		NSUInteger kUnit = 1;
		NSUInteger lUnit = 1;
		NSInteger omegaSign = 0;
		
		/************************************************************************************************/
		/*		Define the problem dimensions															*/
		/************************************************************************************************/
		
		GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:32 domainMin: -L_domain/2 length:L_domain];
		xDim.name = @"x";
		GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:32 domainMin: -L_domain/2  length:L_domain];
		yDim.name = @"y";
		GLMutableDimension *tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
		tDim.name = @"time";
		
		// Variables are always tied to a particular equation---so we create an equation object first.
		GLEquation *equation = [[GLEquation alloc] init];
		
		NSArray *spatialDimensions = @[xDim, yDim];
		GLFunction *x = [GLFunction functionOfRealTypeFromDimension: xDim withDimensions: spatialDimensions forEquation: equation];
//		GLFunction *y = [GLFunction functionOfRealTypeFromDimension: yDim withDimensions: spatialDimensions forEquation: equation];
		
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
		GLFunction *U_mag = [[[omega multiply: @(1/f0)] pow: -1.5] multiply: @(U_max)];
		
		for (NSUInteger i=0; i<kDim.nPoints; i++) {
			for (NSUInteger j=1; j<lDim.nPoints-1; j++) {
				U_mag.pointerValue[ (kDim.nPoints/2)*lDim.nPoints + j] = U_mag.pointerValue[ (kDim.nPoints/2)*lDim.nPoints + j]/2;
			}
		}
		
		// Initial randomized phase
		GLFunction *phi0_plus = [GLFunction functionWithRandomValuesBetween: -M_PI and: M_PI withDimensions: spectralDimensions forEquation: equation];
		GLFunction *phi0_minus = [GLFunction functionWithRandomValuesBetween: -M_PI and: M_PI withDimensions: spectralDimensions forEquation: equation];
		if (shouldUnitTest) {
//			[phi0_plus zero];
//			[phi0_minus zero];
			NSLog(@"omega/k = %f m/s", omega.pointerValue[kUnit*lDim.nPoints+lUnit]/kH.pointerValue[kUnit*lDim.nPoints+lUnit]);
		}
				
		for (NSUInteger i=1; i<kDim.nPoints/2; i++) {
			// Fix the l=0, k=1..N/2-1 components so the phase is the conjugate
			alpha.pointerValue[(kDim.nPoints-i)*lDim.nPoints] = -alpha.pointerValue[i*lDim.nPoints];
			phi0_plus.pointerValue[(kDim.nPoints-i)*lDim.nPoints] = -phi0_plus.pointerValue[i*lDim.nPoints];
			phi0_minus.pointerValue[(kDim.nPoints-i)*lDim.nPoints] = -phi0_minus.pointerValue[i*lDim.nPoints];
			
			// Fix the l=N/2, k=1..N/2-1 components so the phase is the conjugate
			alpha.pointerValue[(kDim.nPoints-i)*lDim.nPoints + lDim.nPoints-1] = -alpha.pointerValue[i*lDim.nPoints + lDim.nPoints-1];
			phi0_plus.pointerValue[(kDim.nPoints-i)*lDim.nPoints + lDim.nPoints-1] = -phi0_plus.pointerValue[i*lDim.nPoints+ lDim.nPoints-1];
			phi0_minus.pointerValue[(kDim.nPoints-i)*lDim.nPoints+ lDim.nPoints-1] = -phi0_minus.pointerValue[i*lDim.nPoints+ lDim.nPoints-1];
		}
		
		// This computes how they map to u,v
		
		GLFunction *fomega = [omega scalarDivide: f0];
		
		GLFunction *u_phase_plus = [[alpha cos] minus: [[[alpha sin] multiply: fomega] swapComplex]];
		GLFunction *v_phase_plus = [[alpha sin] plus: [[[alpha cos] multiply: fomega] swapComplex]];
		
		GLFunction *u_phase_minus = [[alpha cos] plus: [[[alpha sin] multiply: fomega] swapComplex]];
		GLFunction *v_phase_minus = [[alpha sin] minus: [[[alpha cos] multiply: fomega] swapComplex]];
		
		if (shouldUnitTest) {
			if (omegaSign < 0) {
				// We only want the negative rotating part
				[U_mag zero];
				
				// So zero the positive part
				u_phase_plus = [U_mag multiply: u_phase_plus];
				v_phase_plus = [U_mag multiply: v_phase_plus];
				
				// And set the negative part
				U_mag = [U_mag setValue: U_max atIndices: [NSString stringWithFormat:@"%lu,%lu", kUnit, lUnit]];
				
				u_phase_minus = [U_mag multiply: u_phase_minus];
				v_phase_minus = [U_mag multiply: v_phase_minus];
			} else if (omegaSign > 0) {
				// We only want the negative rotating part
				[U_mag zero];
				
				// Zero the negative part
				u_phase_minus = [U_mag multiply: u_phase_minus];
				v_phase_minus = [U_mag multiply: v_phase_minus];
				
				// And set the positive part
				U_mag = [U_mag setValue: U_max atIndices: [NSString stringWithFormat:@"%lu,%lu", kUnit, lUnit]];
				
				u_phase_plus = [U_mag multiply: u_phase_plus];
				v_phase_plus = [U_mag multiply: v_phase_plus];
			} else {
				// We only want the negative rotating part
				[U_mag zero];
				U_mag = [U_mag setValue: U_max atIndices: [NSString stringWithFormat:@"%lu,%lu", kUnit, lUnit]];
				
				// set the negative part
				u_phase_minus = [U_mag multiply: u_phase_minus];
				v_phase_minus = [U_mag multiply: v_phase_minus];
				
				// And set the positive part
				u_phase_plus = [U_mag multiply: u_phase_plus];
				v_phase_plus = [U_mag multiply: v_phase_plus];
			}
            
        } else {
			u_phase_plus = [U_mag multiply: u_phase_plus];
			v_phase_plus = [U_mag multiply: v_phase_plus];
			
			u_phase_minus = [U_mag multiply: u_phase_minus];
			v_phase_minus = [U_mag multiply: v_phase_minus];
		}
		
		
		
		// Fix the l=0, k=1..N/2-1 components
//		for (NSUInteger i=1; i<kDim.nPoints/2; i++) {
//			omega.pointerValue[(kDim.nPoints-i)*lDim.nPoints] = -omega.pointerValue[i*lDim.nPoints];
//			
//			u_phase.pointerValue[(kDim.nPoints-i)*lDim.nPoints] = u_phase.pointerValue[i*lDim.nPoints];
//			u_phase.pointerValue[u_phase.nDataPoints +(kDim.nPoints-i)*lDim.nPoints] = -u_phase.pointerValue[u_phase.nDataPoints + i*lDim.nPoints];
//			
//			v_phase.pointerValue[(kDim.nPoints-i)*lDim.nPoints] = v_phase.pointerValue[i*lDim.nPoints];
//			v_phase.pointerValue[u_phase.nDataPoints +(kDim.nPoints-i)*lDim.nPoints] = -v_phase.pointerValue[u_phase.nDataPoints + i*lDim.nPoints];
//		}
		//		GLFunction *tmp = [u_phase multiply: time_phase];
		
		//		for (NSUInteger i=0; i<u_phase.nDataElements; i++) {
		//			if ( i%lDim.nPoints == 0) printf("\n%0lu: ", i/lDim.nPoints);
		//			printf("%.1g ", u_phase.pointerValue[i]);
		//		}
		
		/************************************************************************************************/
		/*		Now run the model																		*/
		/************************************************************************************************/
		
		
		 NSArray * (^timeToUV) (GLScalar *) = ^( GLScalar *t ) {
			 GLFunction *phi_plus = [phi0_plus plus: [omega multiply: t]];
			 GLFunction *phi_minus = [phi0_minus minus: [omega multiply: t]];
			 
			 GLFunction *time_phase_plus = [[phi_plus swapComplex] exponentiate];
			 GLFunction *time_phase_minus = [[phi_minus swapComplex] exponentiate];
			 
			 GLFunction *u = [[[u_phase_plus multiply: time_phase_plus] plus: [u_phase_minus multiply: time_phase_minus]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis)]];
			 GLFunction *v = [[[v_phase_plus multiply: time_phase_plus] plus: [v_phase_minus multiply: time_phase_minus]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis)]];
			 
			 return @[u,v];
		 };
		
		GLScalar *t = [GLScalar scalarWithValue: 0.0*2*M_PI/f0 forEquation: equation];
		NSArray *uv = timeToUV(t);
		GLFunction *u = uv[0];
		GLFunction *v = uv[1];
        GLFunction *speed = [[[u times: u] plus: [v times: v]] sqrt];
		GLFloat maxSpeed = [speed maxNow];
		NSLog(@"U_max: %f, maxSpeed: %f", U_max, maxSpeed);
		
        /************************************************************************************************/
		/*		Let's also plop a float at each grid point.                                             */
		/************************************************************************************************/
        
        GLDimension *xFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:1 domainMin: 0 length:L_domain];
		xFloatDim.name = @"x-float";
		GLDimension *yFloatDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:1 domainMin: 0  length:L_domain];
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
        
        CGFloat cfl = 0.25;
        GLFloat timeStep = cfl * xDim.sampleInterval / maxSpeed;
        GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: @[xPosition, yPosition] stepSize: timeStep fFromTY:^(GLScalar *time, NSArray *yNew) {
            
			NSArray *uv = timeToUV(time);
			GLFunction *u2 = uv[0];
			GLFunction *v2 = uv[1];
        
			GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: @[u2, v2] secondOperand: yNew];
			
			return interp.result;
		}];
		
//		NSLog(integrator.graphvisDescription);
        
        GLFloat maxTime = maxInertialPeriods*2*M_PI/f0;
        GLFloat sampleTime = sampleTimeInMinutes*60;
        for (GLFloat time = sampleTime; time < maxTime; time += sampleTime)
        {
            @autoreleasepool {
                NSLog(@"Logging day: %f, step size: %f.", (time/86400), integrator.lastStepSize/86400);
                
                NSArray *yout = [integrator stepForwardToTime: time];
                
                t = [GLScalar scalarWithValue: time forEquation: equation];
                NSArray *uv = timeToUV(t);
				GLFunction *u = uv[0];
				GLFunction *v = uv[1];
                
//				GLFunction *xpos = yout[0];
//				GLFunction *ypos = yout[1];
//				NSLog(@"(x,y)=(%f,%f)", xpos.pointerValue[0], ypos.pointerValue[0]);
				
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

