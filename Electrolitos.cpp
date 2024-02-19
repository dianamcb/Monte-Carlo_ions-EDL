//Para este caso, la capa de Stern no depende de la reducción del diámetro de los iones.
//José Guadalupe Ibarra Armenta
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define CLASES 500
#define MAXPART 2100// partículas máximas a manejar
#define pi 3.14159265358979323846
#define kb 1.38062e-23
#define esc 1e-9//más definiciones, para cálculo de la energía 
#define epce 8.85419e-12
#define qe 1.60219e-19
#define na 6.02214e23

void lee_entrada(void);
void imprime_entrada(void);
double alea(void);
void imprime_gdr(void);
void gdr(void);
void asigna_posiciones(void);
void calculos_iniciales(void);
void salida3d(void);
void selecciona_ion(void);
void mueve_ion(void);
void condiciones_periodicas(void);
int signo(float a);
void metropolis(void);
float variacion_energia(void);
float distancia(int a, int b);
int esferas_duras(int a);
void inmovilidad(void);
int esferas_duras_2 (int a);
void energia_total(void);
void inicializando(void);

struct io{//estructura
	float x, y, z;
	int carga, inmovilidad;
}ion [MAXPART];

int ni, n, np, nn, nz=0;
int p, pasos, actu, terma, vp, vn, inmov;
float temp, concentracion, epsi, sigma;
float diam, frac_diam, etotal;//nueva definición: frac_diam//***cambiamos el nombre de la variable para mayor claridad**** Chema
float w,l;
float dl,dw, densp, denspp;
float dx, dy, dz;
int aceptacion=0, rechazo=0, rechazo_traslape=0, rechazo_metropolis=0, rechazo_pared=0, conteo_inmovilidad=0;//nueva definicion: inmovilidad Dy s
long int gdxp[CLASES+1]={0}, gdyp[CLASES+1]={0}, gdzp[CLASES+1]={0}, gdxn[CLASES+1]={0}, gdyn[CLASES+1]={0}, gdzn[CLASES+1]={0};//positivos y negativos por separado
char datos[35];
FILE *dat;
//INICIO DEL PROGRAMA PRINCIPAL***************************************************************************************************************************
main(){
	float diatemp;
	srand((unsigned)time(NULL));/*sembrando la semilla*/
	system ("mkdir temp");
	lee_entrada();
	imprime_entrada();
	//for (diatemp=0.4; diatemp<=1; diatemp+=0.1){
		frac_diam=1;
		calculos_iniciales();
		inicializando();
		printf("\ndiam: %0.2f<------------------------------------------------------------------\n",frac_diam);
		sprintf(datos, "cd temp & mkdir diam_%0.3f", frac_diam);
		system(datos);
		sprintf(datos, "temp/diam_%0.3f/energia.dat", frac_diam);
		dat=fopen(datos, "w");
		fclose(dat);	
		asigna_posiciones();
		for(p=1; p<=pasos; p++){/*INICIO(Principal)*/
			rechazo=0;
			selecciona_ion();
			mueve_ion();
			condiciones_periodicas();
			if(rechazo==0) metropolis();		
			if(ion[ni].inmovilidad>=inmov){//La func. inmovilidad mueve al ión a otra posición al azar.	
				inmovilidad();
				conteo_inmovilidad++;
			}
			if(p>terma-actu) gdr();
			//Se muestra el promedio de inmovilidad en la consola junto con el promedio de los rechazos. D y S
			if(p%10000==0) printf("\rPaso: %i Acep=%0.3f r_met: %0.3f r_tras: %0.3f, r_pared: %0.3f, inmovi: %0.3f, UNO: %0.3f n: %i", p, aceptacion*1.0/p, rechazo_metropolis*1.0/p, rechazo_traslape*1.0/p, rechazo_pared*1.0/p, conteo_inmovilidad*1.0/p, (rechazo_pared*1.0+aceptacion*1.0+rechazo_traslape*1.0+rechazo_metropolis*1.0)/p, n);
			if((p%actu==0)&&(p>=terma)){
				//energia_total();
				imprime_gdr();
				salida3d();
			}
		}/*FINAL(Principal)*/
	//}
	getchar();
	return(0);
}
//FINAL DEL PROGRAMA PRINCIPAL*****************************************************************************************************************************
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double alea(void){//genera números aleatorios entre 0 y 1
return((double)rand()/RAND_MAX); //32767 dependiente de la libreria
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void inicializando(void){
	int i;
	aceptacion=rechazo=rechazo_traslape=rechazo_metropolis=rechazo_pared=0;//nuevas definiciones
	for(i=1; i<=CLASES+1; i++){
		gdxp[i]=gdyp[i]=gdzp[i]=gdxn[i]=gdyn[i]=gdzn[i]=0;//positivos y negativos por separado
	}
	conteo_inmovilidad=0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void gdr(void){//modificaciones
	int i;
	for(i=1; i<=n; i++){
		if(ion[i].carga==vp){
			gdxp[int(ion[i].x/dw)+1]++;	
			gdyp[int(ion[i].y/dw)+1]++;	
			gdzp[int((ion[i].z-frac_diam/2)/dl)+1]++;//
		}else{
			gdxn[int(ion[i].x/dw)+1]++;	
			gdyn[int(ion[i].y/dw)+1]++;	
			gdzn[int((ion[i].z-frac_diam/2)/dl)+1]++;//	
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void imprime_gdr(void){//modificaciones
	int i, pr=p-terma+actu;
	float dens_carga=-sigma, area=pow(w*diam*esc, 2), potencial=0;//area es el área de la cara superior de las subplacas formadas por fraccionar el volumen neto por CLASES.
	FILE *dat;
	char nombre[50];
	sprintf(nombre,"temp/diam_%0.3f/gdr_%i.dat", frac_diam, p/actu);
	//sprintf(nombre,"gdr.dat");
	dat=fopen(nombre, "w");
		for(i=1; i<=CLASES; i++){
			potencial+=(frac_diam/2 - ((i-0.5)*dl+frac_diam/2))*(esc*diam)*(gdzp[i]*vp-gdzn[i])*qe*dl*(esc*diam)/(epce*epsi*pr*w*w*dl*pow(esc*diam,3));
			dens_carga+=1e2*(gdzp[i]*vp-gdzn[i])*qe/(pr*area);
			fprintf(dat,"%f	%f	%f	%f	%f	%f	%f	%f	%f	%f\n", (i-0.5)*dw, gdxp[i]/(1.0*pr*denspp), gdyp[i]/(1.0*pr*denspp), gdxn[i]/(1.0*pr*densp), gdyn[i]/(1.0*pr*densp), (i-0.5)*dl+frac_diam/2, gdzp[i]*vp/(1.0*pr*densp), gdzn[i]/(1.0*pr*densp), dens_carga, potencial);//====vp===distribucion en las tres componentes de r
		}
	fclose(dat);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void asigna_posiciones(void){
int i;
float xi, fxi;//declaramos i y novimos las otras declaraciones
	for(i=1; i<=n; i++){
		if (i<=np) ion[i].carga=vp;
		else ion[i].carga=vn;
		//do{
			xi=alea()*w;
			fxi=exp(-pow(xi-w/2,2)/2)/sqrt(2*pi);
		//}while(alea()>fxi);
		ion[i].x=xi;
		//do{
			xi=alea()*w;
			fxi=exp(-pow(xi-w/2,2)/2)/sqrt(2*pi); 
		//}while(alea()>fxi);
		ion[i].y=xi;
		//do{
			xi=alea()*(l-frac_diam)+frac_diam/2;
			//fxi=exp(-pow(xi-l/2,2)/2)/sqrt(2*pi);
		//}while(alea()>fxi);
		ion[i].z=xi;
		if(esferas_duras(i)==1) i--;
	}	
	/*ion[n+1].x=10;
	ion[n+1].y=5;
	ion[n+1].z=10;
	ion[n+1].carga=200;
	ion[n+2].x=10;
	ion[n+2].y=15;
	ion[n+2].z=10;
	ion[n+2].carga=-200;*/
	salida3d();
	p=n;//truco;
	gdr();//importante llamada
	imprime_gdr();//salida 0
	/*gdxp[CLASES+1]={0};//reinicializar estadisticas
	gdyp[CLASES+1]={0};
	gdzp[CLASES+1]={0};
	gdxn[CLASES+1]={0};
	gdyn[CLASES+1]={0};
	gdzn[CLASES+1]={0};*/
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void salida3d(void){
	int m;
	FILE *dat;
	char nombre[30];
	sprintf(nombre, "temp/diam_%0.3f/salida3d%i.dat", frac_diam, p/actu);
	dat=fopen(nombre, "w");
	fprintf(dat, "Execute[{");
	for (m=1; m<=n;m++){
		if (frac_diam==0){
			if(ion[m].carga==vp) fprintf(dat,"\"S%i:(x-%f)^2+(y-%f)^2+(z-%f)^2=%f\",\"SetColor[S%i,Blue]\"", m, ion[m].x, ion[m].y, ion[m].z, pow(0.1/2,2), m);
			else if (ion[m].carga==vn) fprintf(dat,"\"S%i:(x-%f)^2+(y-%f)^2+(z-%f)^2=%f\",\"SetColor[S%i,Red]\"", m, ion[m].x, ion[m].y, ion[m].z, pow(0.1/2,2), m);
		}else{
			if(ion[m].carga==vp) fprintf(dat,"\"S%i:(x-%f)^2+(y-%f)^2+(z-%f)^2=%f\",\"SetColor[S%i,Blue]\"", m, ion[m].x, ion[m].y, ion[m].z, pow(frac_diam/2,2), m);
			else if (ion[m].carga==vn) fprintf(dat,"\"S%i:(x-%f)^2+(y-%f)^2+(z-%f)^2=%f\",\"SetColor[S%i,Red]\"", m, ion[m].x, ion[m].y, ion[m].z, pow(frac_diam/2,2), m);
		}
		//else   fprintf(dat,"\"S%i:(x-%f)^2+(y-%f)^2+(z-%f)^2=%f\",\"SetColor[S%i,Yellow]\"", m, ion[m].x, ion[m].y, ion[m].z, pow(diam/2,2), m);
		if(m<n)fprintf(dat, ",");
	}
	fprintf(dat, "}]");
	fclose(dat);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void selecciona_ion(void){
	do{
		ni=int(alea()*n)+1;
	}while(ni>n);
//printf("seleccionar\n");
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void mueve_ion(void){
	ion[0]=ion[ni];
	ion[ni].x+=(alea()-0.5)*dx;
	ion[ni].y+=(alea()-0.5)*dy;
	ion[ni].z+=(alea()-0.5)*dz;
	/*ion[ni].x=(alea())*w;
	ion[ni].y=(alea())*w;
	ion[ni].z=(alea())*l;*/
//	printf("mover\n");
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void condiciones_periodicas(void){//check backwards
	if((ion[ni].x<0) || (ion[ni].x>=w)){
		ion[ni].x+=w*signo(-ion[ni].x);
//		printf("condiciones.x\n");
	}
	if((ion[ni].y<0) || (ion[ni].y>=w)){
		ion[ni].y+=w*signo(-ion[ni].y);
//		printf("condiciones.y\n");
	}
	/*if(sigma==0){
		if((ion[ni].z<0.0) || (ion[ni].z>=l))   ion[ni].z-=l*signo(ion[ni].z);
	}else{*/
		if((ion[ni].z<frac_diam/2) || (ion[ni].z>=(l-frac_diam/2))){
	//		ion[ni].z+=l*signo(-ion[ni].z);
			rechazo_pared++;
			ion[ni].inmovilidad++;
			rechazo=1;
			ion[ni]=ion[0];//Se aplican condiciones periódicas en "x" y "y", pero no en z.
	//		printf("condiciones.z\n");
		}
	//}
//*************************************************NUEVO*************************************************************
	//if((ion[ni].z<(frac_diam/2.0)) || (ion[ni].z>=l-(frac_diam/2.0)){
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int signo(float a){
	if(a>=0) return(1);
	else if (a<0) return(-1);
	//else return (0);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void metropolis(void){
//	printf("metropolis.inicio\n");
//	getchar();
	double argexp;
		argexp=(-1)*variacion_energia()/(kb*temp);
	if(rechazo==0){
		if (exp(argexp)>alea()){
			aceptacion++;
			ion[ni].inmovilidad=0;
			//printf"
		}else{
			ion[ni]=ion[0];
			rechazo_metropolis++;
			ion[ni].inmovilidad++;
		}
	}else{//el programa solo entra aca si hubo traslape de esferas
		ion[ni]=ion[0];
		rechazo_traslape++;
		ion[ni].inmovilidad++;
	}
//	printf("metropolis\n");
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float variacion_energia(void){
	int i;
	float ecoul_i=0, ecoul_f=0, eplano_i=0, eplano_f=0, den, cte=qe*qe/(4*pi*epce*epsi), lr_i=0, lr_f=0, cte_lr=-qe*qe/((pi*epce*epsi)*pow(w*diam*esc,2)), z_ij, r_1, r_2;//aquí introduje la permitividad dielectrica relativa del agua a 300 K para suavizar las interacciones
	eplano_i=sigma*1e-2*ion[0].carga*qe*(ion[0].z*diam*esc)/(2*epce*epsi);
	eplano_f=sigma*1e-2*ion[ni].carga*qe*(ion[ni].z*diam*esc)/(2*epce*epsi);
	for (i=1; i<=n;i++){
		if(i!=ni){
			if (distancia(ni,i)>=frac_diam){
				ecoul_i+=cte*ion[0].carga*ion[i].carga/(distancia(0,i)*diam*esc);
				ecoul_f+=cte*ion[ni].carga*ion[i].carga/(distancia(ni,i)*diam*esc);
				z_ij=fabs(ion[0].z-ion[i].z);
				r_1=sqrt(0.5+pow(z_ij/w,2));
				r_2=sqrt(0.25+pow(z_ij/w,2));
				lr_i+=ion[0].carga*ion[i].carga*cte_lr*((w*diam*esc)*log((0.5+r_1)/r_2)+(z_ij*diam*esc)*atan(4*z_ij*r_1/w));//Hay que revisar los factores de escala
				z_ij=fabs(ion[ni].z-ion[i].z);
				r_1=sqrt(0.5+pow(z_ij/w,2));
				r_2=sqrt(0.25+pow(z_ij/w,2));
				lr_f+=ion[0].carga*ion[i].carga*cte_lr*((w*diam*esc)*log((0.5+r_1)/r_2)+(z_ij*diam*esc)*atan(4*z_ij*r_1/w));//Hay que revisar los factores de escala
	    	}else {//en esta parte se revisa si los iones se traslapan y de ser asi se rechaza el movimiento se corta el calculo y la variable rechazo cambia de valor para que no se efectue el algoritmo de Metropolis
				rechazo=1;
//				printf("traslape\n");
				return (den);
			}
		}
	}
	den=ecoul_f-ecoul_i+eplano_f-eplano_i+lr_f-lr_i;
//	printf("\necoul_i: %e, ecoul_f: %e \neplano_i: %e, eplano_f: %e \nlr_i: %e, lr_f: %e \nden: %e \n", ecoul_i, ecoul_f, eplano_i, eplano_f, lr_i, lr_f, den);
//	printf("Den=%1.6e\n",den);
//	printf("variacion\n");
	return(den);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void energia_total(void){
	int i, j;
	float z_ij, r_1, r_2;
	float eplano=0, ecoul=0, lr=0, cte=qe*qe/(4*pi*epce*epsi), cte_lr=-qe*qe/((pi*epce*epsi)*pow(w*diam*esc,2));//aquí se introdujo la permitividad dielectrica relativa del agua a 300 K para suavizar las interacciones.
	char nombre[30];
	for(i=1; i<=n; i++){
		eplano+=sigma*1e-2*ion[i].carga*qe*(ion[i].z*diam*esc)/(2*epce*epsi);
		for(j=1; j<=n; j++){
			if(j>i){
				ecoul+=cte*ion[i].carga*ion[j].carga/(distancia(i,j)*diam*esc);
				z_ij=fabs(ion[i].z-ion[j].z);
				r_1=sqrt(0.5+pow(z_ij/w,2));
				r_2=sqrt(0.25+pow(z_ij/w,2));
				lr+=ion[i].carga*ion[j].carga*cte_lr*((w*diam*esc)*log((0.5+r_1)/r_2)+(z_ij*diam*esc)*atan(4*z_ij*r_1/w));//Hay que revisar los factores de escala.
			}
		}		
	}
	etotal=ecoul+eplano+lr;
	/*Se registrará la energía por cada vez que se imprima un archivo temp/gdr_%i dentro del archivo energia.dat*/
	sprintf(nombre, "temp/diam_%0.3f/energia.dat", frac_diam);
	dat=fopen(nombre, "a");
		fprintf(dat, "%e\n", etotal);
	fclose(dat);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float distancia(int a, int b){
	float dist, xx, yy, zz;
	xx=ion[b].x;
	yy=ion[b].y;
	zz=ion[b].z;
	if (fabs(ion[a].x-xx)>w/2) xx+=w*signo(ion[a].x-xx);//check backwards
	if (fabs(ion[a].y-yy)>w/2) yy+=w*signo(ion[a].y-yy);
	/*if(sigma==0){
		if(fabs(ion[a].z-zz)>l/2)zz+=l*signo(ion[a].z-zz);
	}*/
	dist=sqrt(pow(ion[a].x-xx,2)+pow(ion[a].y-yy,2)+pow(ion[a].z-zz,2));//poner math.h
	return (dist);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int esferas_duras(int a){
	int i, res=1;
	for(i=1; i<=a-1; i++){
		if(distancia(i,a)<frac_diam) return(res);
	}
	res=0;
	return(res);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int esferas_duras_2(int a){
	int i, res=1;
	for(i=1; i<=n; i++) {
		if(distancia(i,a)<frac_diam && i!=a) return(res);
	}
	res=0;
	return(res);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lee_entrada(void){
	dat=fopen("entrada.dat","r");
	fscanf(dat,"pasos: %i", &pasos);
	fscanf(dat,"\ntermalizacion: %i", &terma);
	fscanf(dat,"\nactualizacion: %i", &actu);
	fscanf(dat,"\ninmovilidad: %i", &inmov);
	fscanf(dat,"\nvalencias(v+,v-): %i, %i", &vp, &vn);
	fscanf(dat,"\nconcentracion(M/lt): %f", &concentracion);
	fscanf(dat,"\ntemperatura(k): %f", &temp);
	fscanf(dat,"\npermitividad dielectrica relativa: %f", &epsi);
	fscanf(dat,"\ndensidad de carga(uC/cm^2): %f", &sigma);
	fscanf(dat,"\ndiametro(nm): %f", &diam);
	fscanf(dat,"\ndiametro reducido(fraccion del diametro): %f", &frac_diam);
	fscanf(dat,"\n(w,l)(diametros): %f, %f", &w, &l);
	fscanf(dat,"\ndx, dy, dz: %f, %f, %f", &dx, &dy, &dz);
	fclose(dat);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void imprime_entrada(void){
	printf("pasos: %i", pasos);
	printf("\ntermalizacion: %i", terma);
	printf("\nactualizacion: %i", actu);
	printf("\ninmovilidad: %i", inmov);
	printf("\nvalencias(v+,v-): %i, %i", vp, vn);
	printf("\nconcentracion(M/lt): %f", concentracion);
	printf("\ntemperatura(k): %f", temp);
	printf("\npermitividad dielectrica relativa: %f", epsi);
	printf("\ndensidad de carga(uC/cm^2): %f", sigma);//Aunque la densidad sea negativa, debe anotarse como positiva.
	printf("\ndiametro(nm): %f", diam);
	printf("\ndiametro reducido(fraccion del diametro): %f", frac_diam);
	printf("\n(w,l)(diametros): %f %f", w, l);
	printf("\ndx,dy,dz: %f, %f, %f", dx, dy, dz);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculos_iniciales(void){
	float wr;
	if (sigma!=0){
		nz=int((sigma*1e-2*pow(w*diam*esc,2))/(qe*vp)); //Se calcula el número de cargas en el plano z.
		wr=sqrt(qe*nz*vp/(sigma*1e-2*pow(diam*esc,2)));//Definimos a wr para obtener un valor nz no decimal (si no contamos el comando int())
		w=wr;
	}
	dl=(l-frac_diam)/CLASES;
	dw=w/CLASES;
	nn=vp*int(concentracion*1000*na*pow(w*diam*esc, 2)*(l-frac_diam)*diam*esc);//vp*
	np=(nn/vp)+nz;
	n=nn+np;
	denspp=1.0*np/CLASES;
	densp=nn*1.0/CLASES;
	printf("\nwr: %f", wr);
	printf("\ndl: %f, dw: %f", dl, dw);
	printf("\nnn, np, nz, n: %i, %i, %i, %i\n", nn, np, nz, n);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void inmovilidad(void){
	do{
		ion[ni].x=alea()*w;
		ion[ni].y=alea()*w;
		//ion[ni].z=alea()*(l-1)+0.5;
	}while(esferas_duras_2(ni)==1 || (ion[ni].x==w && ion[ni].y==w));
	ion[ni].inmovilidad=0;
}
/////////////////////////////////////////////////////////////////////////////THE END///////////////////////////////////////////////////////////////////////	
