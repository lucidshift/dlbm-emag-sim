void Curl(){

	double B2, H2, B2mag, H2mag;
	double muTemp;

    for (int x=1;x<XDIM-1;x++){
      	for (int y=1;y<YDIM-1;y++){
			for(int z=1; z<(ZDIM-1); z++){

				muTemp = u;
				
				//B field is calculated and redimensionalized here.
      			B[x][y][z][0] = (4 * pi * ((rho[2][x][y+1][z] - rho[2][x][y-1][z]) - (rho[1][x][y][z+1] - rho[1][x][y][z-1]))) / 2;
      			B[x][y][z][1] = (4 * pi * ((rho[0][x][y][z+1] - rho[0][x][y][z-1]) - (rho[2][x+1][y][z] - rho[2][x-1][y][z]))) / 2;
      			B[x][y][z][2] = (4 * pi * ((rho[1][x+1][y][z] - rho[1][x-1][y][z]) - (rho[0][x][y+1][z] - rho[0][x][y-1][z]))) / 2;

				B2 = sqrt(((B[x][y][z][0]*B[x][y][z][0]) + (B[x][y][z][1]*B[x][y][z][1]) + (B[x][y][z][2]*B[x][y][z][2])));	

				H[x][y][z][0] = B[x][y][z][0] / MPERM;
				H[x][y][z][1] = B[x][y][z][1] / MPERM;
				H[x][y][z][2] = B[x][y][z][2] / MPERM;

				H2 = sqrt((H[x][y][z][0]*H[x][y][z][0]) + (H[x][y][z][1]*H[x][y][z][1]) + (H[x][y][z][2]*H[x][y][z][2]));


				//This is the H field, because the u0 factor was moved to the other side, multiplied into the JInduced current. Divided by two due to how the difference is calculated.
      			Hmag[x][y][z][0] = (4 * pi * ((rhomag[2][x][y+1][z] - rhomag[2][x][y-1][z]) - (rhomag[1][x][y][z+1] - rhomag[1][x][y][z-1]))) / 2;
      			Hmag[x][y][z][1] = (4 * pi * ((rhomag[0][x][y][z+1] - rhomag[0][x][y][z-1]) - (rhomag[2][x+1][y][z] - rhomag[2][x-1][y][z]))) / 2;
      			Hmag[x][y][z][2] = (4 * pi * ((rhomag[1][x+1][y][z] - rhomag[1][x-1][y][z]) - (rhomag[0][x][y+1][z] - rhomag[0][x][y-1][z]))) / 2;

      			H2mag = sqrt(((Hmag[x][y][z][0]*Hmag[x][y][z][0]) + (Hmag[x][y][z][1]*Hmag[x][y][z][1]) + (Hmag[x][y][z][2]*Hmag[x][y][z][2])));

      			Bmag[x][y][z][0] = Hmag[x][y][z][0] / mu[x][y][z];
      			Bmag[x][y][z][1] = Hmag[x][y][z][1] / mu[x][y][z];
      			Bmag[x][y][z][2] = Hmag[x][y][z][2] / mu[x][y][z];

				B2mag = sqrt(((Bmag[x][y][z][0]*Bmag[x][y][z][0]) + (Bmag[x][y][z][1]*Bmag[x][y][z][1]) + (Bmag[x][y][z][2]*Bmag[x][y][z][2])));


				/*
				//Permeability Curve
				if(B2<0.0 && muTemp != u0){		//Don't change the permeability of free space, it doesn't change.
					mu[x][y][z] = u;
				}
				if(B2>0.0 && B2<1.5 && muTemp != u0){
					mu[x][y][z] = (3484.714 * (B2*B2*B2*B2) - 6261.368 * (B2*B2*B2) - 9598.177 * (B2*B2) + 15953.587 * B2 + 3119.105);
				}	
				if(B2>=1.5 && muTemp != u0){
					mu[x][y][z] = 608;
				}	
				*/								

				//H = B/u0 - M
				Htotal[x][y][z][0] = H[x][y][z][0] - Hmag[x][y][z][0];
				Htotal[x][y][z][1] = H[x][y][z][1] - Hmag[x][y][z][1];
				Htotal[x][y][z][2] = H[x][y][z][2] - Hmag[x][y][z][2];


				if(x==27 && y==27 && z==50){
					//printf("|B|=%10.10f\n", B2);
					//printf("Mu=%10.10f\n",muTemp);

					BMiddle = B2;

					if(((iterationCount*2)+1)<sizeX){
						BHCurve[iterationCount*2] = H2;	
						BHCurve[iterationCount*2 + 1] = B2;
						HUrCurve[iterationCount*2] = H2;
						HUrCurve[iterationCount*2 + 1] = muTemp;
						printf("H2=%5.25f,B2=%5.25f\n", H2, B2);
						printf("H2mag=%5.25f,B2mag=%5.25f\n", H2mag, B2mag);
					}					
				}	

      		}
		}
  	}
}

void GetFieldEnergy(){

	double eTotal = 0;

	for(int x=0; x<XDIM; x++){
		for(int y=0; y<YDIM; y++){
			for(int z=0; z<ZDIM; z++){
			
							// 		sum(B^2 * dv * 1/(2u0))															//dv is 1mm, hence the factor of 10^-9 for the B field
				eTotal += (((B[x][y][z][0] * B[x][y][z][0]) + (B[x][y][z][1]*B[x][y][z][1]) + (B[x][y][z][2]*B[x][y][z][2])) * 0.000000001) / (2 * MPERM);


			}
		}
	}
	Energy[iterationCount] = eTotal;
	printf("Total Field Energy = %10.15f\n", eTotal);
}

void CalcOrbitalCurrent(){

	//double muFactor; 
	double dxMu, dyMu, dzMu;
	//double badMu = - 1./(u0 * u0 * MPERM * MPERM);

	double dxRhoX, dxRhoY, dxRhoZ;
	double dyRhoX, dyRhoY, dyRhoZ;
	double dzRhoX, dzRhoY, dzRhoZ;		

	for(int x=1; x<XDIM-1; x++){
		for(int y=1; y<YDIM-1; y++){
			for(int z=1; z<ZDIM-1; z++){

				//muFactor = (mu[x][y][z] * mu[x][y][z] * MPERM * MPERM) / (1000 * 0.148809);

				/*
				if(iterationCount > 50 && x == 26){
					printf("muFactor=%10.10f\n",muFactor);
				}				
				*/
			
				dxMu = (mu[x-1][y][z] - mu[x+1][y][z]) / 2;
				dyMu = (mu[x][y-1][z] - mu[x][y+1][z]) / 2;
				dzMu = (mu[x][y][z-1] - mu[x][y][z+1]) / 2;

				/*
				if(dxMu > 0 && muFactor == badMu){
					muFactor = -1;//- 1./(mu[x-1][y][z] * mu[x-1][y][z]);
				}

				if(dyMu > 0 && muFactor == badMu){
					muFactor = -1;//- 1./(mu[x][y-1][z] * mu[x][y-1][z]);
				}

				if(dzMu > 0 && muFactor == badMu){
					muFactor = -1;//- 1./(mu[x][y][z-1] * mu[x][y][z-1]);
				}
				*/

				dxRhoX = (rho[0][x-1][y][z] - rho[0][x+1][y][z]) / 2;
				dxRhoY = (rho[1][x-1][y][z] - rho[1][x+1][y][z]) / 2;
				dxRhoZ = (rho[2][x-1][y][z] - rho[2][x+1][y][z]) / 2;

				dyRhoX = (rho[0][x][y-1][z] - rho[0][x][y+1][z]) / 2;
				dyRhoY = (rho[1][x][y-1][z] - rho[1][x][y+1][z]) / 2;
				dyRhoZ = (rho[2][x][y-1][z] - rho[2][x][y+1][z]) / 2;

				dzRhoX = (rho[0][x][y][z-1] - rho[0][x][y][z+1]) / 2;
				dzRhoY = (rho[1][x][y][z-1] - rho[1][x][y][z+1]) / 2;
				dzRhoZ = (rho[2][x][y][z-1] - rho[2][x][y][z+1]) / 2;	

				//Modified to test uJ hypothesis
				JInduced[0][x][y][z] = -(mu[x][y][z]/u0) * ((dxMu * dxRhoX) + (dyMu * dxRhoY) + (dzMu * dxRhoZ));// * muFactor;
				JInduced[1][x][y][z] = -(mu[x][y][z]/u0) * ((dxMu * dyRhoX) + (dyMu * dyRhoY) + (dzMu * dyRhoZ));// * muFactor;
				JInduced[2][x][y][z] = -(mu[x][y][z]/u0) * ((dxMu * dzRhoX) + (dyMu * dzRhoY) + (dzMu * dzRhoZ));// * muFactor;	

				/*
				if(iterationCount > 50 && x == 26){
					printf("(%10.10f, %10.10f, %10.10f)\n", JInduced[0][x][y][z], JInduced[1][x][y][z], JInduced[2][x][y][z]);
				}
				*/
						
			}
		}
	}
}