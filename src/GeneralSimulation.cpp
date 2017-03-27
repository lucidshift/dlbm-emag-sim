void Initialize(){

	//Sets up u0
	ConfigureFreeSpace();

	//Coil Setup - This section needs work.
	for(int z=57; z<(123); z=z+3){		//Fixed number of coils!!
		int r=11;
		//for(int r=11; r<25; r=r+3){

			coilMidpoint = ((123 - 57) / 2) + 57;	//Change this when you change the coil dimensions!
			printf("(%d,%d)\n",z,r);
			GenerateSourcePath(40, 40, z, r, 3);
		//}
	}

	//Current setup
	for(int x=0; x<(XDIM); x++){
		for(int y=0; y<(YDIM); y++){
			for(int z=0; z<(ZDIM); z++){
												
				Source[0][x][y][z] = JMagnitude * J[x][y][z].x * MPERM;
				Source[1][x][y][z] = JMagnitude * J[x][y][z].y * MPERM;
				Source[2][x][y][z] = JMagnitude * J[x][y][z].z * MPERM;
			}
		}
	}	

	if(netForceSim == 0){	
		//Buffer Initializations
		for(int s=0; s<3; s++){
			for(int x=0; x<(XDIM); x++){
				for(int y=0; y<(YDIM); y++){
					for(int z=0; z<(ZDIM); z++){
						for(int i=0; i<velCount; i++){

							A[s][i][x][y][z] = 0; //velDist[i] + Source[i][x][y];

							leftBuffer[s][i][y][z] = 0; 	//Velocity 1
							rightBuffer[s][i][y][z] = 0;	//Velocity 2
							topBuffer[s][i][x][z] = 0; 		//Velocity 3
							bottomBuffer[s][i][x][z] = 0;	//Velocity 4
							frontBuffer[s][i][x][y] = 0; 	//Velocity 5
							backBuffer[s][i][x][y] = 0;		//Velocity 6					
						}
					}
				}
			}	
		}	

		//Density Initialization
		for(int s=0; s<3; s++){
		    for (int x=0;x<(XDIM);x++){
		      	for (int y=0;y<(YDIM);y++){
					for(int z=0; z<(ZDIM); z++){
				      		rho[s][x][y][z] = 0;
					}
				}
		  	}
		}
	}

		for(int s=0; s<3; s++){
			for(int x=0; x<(XDIM); x++){
				for(int y=0; y<(YDIM); y++){
					for(int z=0; z<(ZDIM); z++){

						rhomag[s][x][y][z] = 0;
						for(int i=0; i<velCount; i++){

							Amag[s][i][x][y][z] = 0; //velDist[i] + Source[i][x][y];

							leftBuffermag[s][i][y][z] = 0; 	//Velocity 1
							rightBuffermag[s][i][y][z] = 0;	//Velocity 2
							topBuffermag[s][i][x][z] = 0; 		//Velocity 3
							bottomBuffermag[s][i][x][z] = 0;	//Velocity 4
							frontBuffermag[s][i][x][y] = 0; 	//Velocity 5
							backBuffermag[s][i][x][y] = 0;		//Velocity 6					
						}
					}
				}
			}	
		}

  	//Graph Initialization
	for(int x=0; x<(XDIM); x++){

		BHCurve[x] = 0;
		HUrCurve[x] = 0;		
	}

	//Reset to prevent initialization loop on false solution convergence.
	for(int i=0; i<20; i++){
		EnergyDelta[i] = 0;
	}
	EDeltaPointer = 0;
	EDeltaMax = 0;

	//Initialize Force arrays
    for (int x=0;x<(XDIM);x++){
      	for (int y=0;y<(YDIM);y++){
			for(int z=0; z<(ZDIM); z++){

				Force[x][y][z][0] = 0;
				Force[x][y][z][1] = 0;
				Force[x][y][z][2] = 0;
			}
		}
  	}

	//System variable setup
	iterationCount = 0;
	cpuCount = sysconf(_SC_NPROCESSORS_ONLN);
	printf("Number of CPU cores availible = %d\n", cpuCount);

}

void GenerateSourcePath(int cx, int cy, int cz, int radius, int thickness){

	int dx;
	int dy;		
	int xPos;
	int yPos;

	//Create the "wire" for current to flow through
	for(double theta=0; theta<(2*pi); theta+=.01){
		for(int r=(radius-(thickness/2)); r<(radius+((thickness+1)/2)); r++){
			for(int z=(cz-1); z<(cz+2); z++){

				dx = cos(theta) * r;
				dy = sin(theta) * r;
				xPos = cx + dx;
				yPos = cy + dy;

				//2D, 90 degree rotation in the x-y plane
				J[xPos][yPos][z].x = -dy;
				J[xPos][yPos][z].y = dx;
				J[xPos][yPos][z].z = 0;
			}				
		}
	}
}

void ConfigurePermeability(double a){

	int cx = 40;
	int cy = 40;
	//int cz = 25;
	int radius = 9;

	int dx;
	int dy;		
	int xPos;
	int yPos;

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
				mu[x][y][z] = u0;
			}
		}
	}

	
	for(double theta=0; theta<(2*pi); theta+=.1){
		for(double r=0; r<(radius); r+=.5){
			for(int z=25; z<ZDIM-25; z++){

				dx = cos(theta) * r;
				dy = sin(theta) * r;
				xPos = cx + dx;
				yPos = cy + dy;

				mu[xPos][yPos][z] = a;

			}				
		}
		//printf("Rad= %5.5f\n",theta);
	}
}

void ConfigureFreeSpace(){

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
				mu[x][y][z] = u0;
			}
		}
	}
}

void ConfigureProjectile(int midpoint, int length, int radius, int permeability){


	int cx = 40;
	int cy = 40;
	//int cz = 25;

	int dx;
	int dy;		
	int xPos;
	int yPos;

	if(midpoint-(length/2) < 0){
		printf("Projectile out of lower simulation bounds.\n");
	}

	if(midpoint+(length/2) >= ZDIM){
		printf("Projectile out of upper simulation bounds.\n");
	}	
	
	for(double theta=0; theta<(2*pi); theta+=.1){
		for(double r=0; r<(radius); r+=.5){
			for(int z=midpoint-(length/2); z<midpoint+(length/2); z++){

				dx = cos(theta) * r;
				dy = sin(theta) * r;
				xPos = cx + dx;
				yPos = cy + dy;

				mu[xPos][yPos][z] = permeability;

			}				
		}
		//printf("Rad= %5.5f\n",theta);
	}


}


void CalcForce(){

	double uTemp = 0;
	double forceTemp[3] = {0,0,0};

	//Potential Energy Calculation
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){

				if(mu[x][y][z] != u0){
					Umag[x][y][z] = ((Hmag[x][y][z][0]*B[x][y][z][0]) + (Hmag[x][y][z][1]*B[x][y][z][1]) + (Hmag[x][y][z][2]*B[x][y][z][2]));
					uTemp += Umag[x][y][z];
				}
			}
		}
	}

	//Force Calculation
	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){

				if(mu[x][y][z] != u0){
					Force[x][y][z][0] = (Umag[x+1][y][z] - Umag[x-1][y][z]) / 2;
					Force[x][y][z][1] = (Umag[x][y+1][z] - Umag[x][y-1][z]) / 2;
					Force[x][y][z][2] = (Umag[x][y][z+1] - Umag[x][y][z-1]) / 2;

					forceTemp[0] += Force[x][y][z][0];
					forceTemp[1] += Force[x][y][z][1];
					forceTemp[2] += Force[x][y][z][2];
				}
			}
		}
	}	

	UmagTotal = uTemp;
	ForceTotal[0] = forceTemp[0];
	ForceTotal[1] = forceTemp[1];
	ForceTotal[2] = forceTemp[2];

	//For multiposition simulations
	ForceVsPosition[GlobalXYZ[2]][0] = forceTemp[0];
	ForceVsPosition[GlobalXYZ[2]][1] = forceTemp[1];
	ForceVsPosition[GlobalXYZ[2]][2] = forceTemp[2];
	ZForceDisplay[GlobalXYZ[2]] = forceTemp[2];

	printf("UmagTotal=%5.15f\n", UmagTotal);
	printf("ForceTotal=(%5.15f, %5.15f, %5.15f)\n", forceTemp[0], forceTemp[1], forceTemp[2]);
}