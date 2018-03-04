//Bibliotecas
//parâmetros leitura analogica
int a0, a1, a2, a3, a4;
long i = 0;
///parâmetros do filtro de kalman
const float L_1[4] = {0.093068029211568, -0.018048815013733, 0.001381271184861, -0.027780007093974};
const float L_2[4] = {0.093068029211568, -0.018048815013733, 0.001381271184861, -0.027780007093974};
float L1[4], L2[4], L3[4], L4[4];
float Ad[4][4];
float C[4];
float B[4], B2[4], B3[4], B4[4];
float FK[4][4], FK2[4][4];
const float Cd[4] = {1, 0, 0, 0};
float xa1[4] = {0, 0, 0, 0}, xa2[4] = {0, 0, 0, 0}, xa3[4] = {0, 0, 0, 0}, xa4[4]  = {0, 0, 0, 0};
String gg1, gg2;
//Parâmetro MAV
long max1 = 0, max2 = 0, max3 = 0, max4 = 0, max5 = 0, max6 = 0, mix1 = 0, mix2 = 0, mix3 = 0, mix4 = 0, mix5 = 0, mix6 = 0;
int WINDOW = 2000, j;

// Covariança do fusor12
float B12[4], xa12[4]  = {0, 0, 0, 0}, P_[4][4], invP_[4][4], Pa[4][4], Pi[4][4];
float  m1[4] = {0, 0, 0, 0}, m2[4] = {0, 0, 0, 0}, M1[4] = {0, 0, 0, 0}, M2[4] = {0, 0, 0, 0}, XA1[4] = {0, 0, 0, 0}, XA2[4] = {0, 0, 0, 0};
//float m1,m2;

// Covariança do fusor34
float B34[4], xa34[4]  = {0, 0, 0, 0}, P_1[4][4],P_2[4][4], invP_1[4][4], Pa1[4][4], Pi1[4][4];
float  m3[4] = {0, 0, 0, 0}, m4[4] = {0, 0, 0, 0}, M3[4] = {0, 0, 0, 0}, M4[4] = {0, 0, 0, 0}, XA3[4] = {0, 0, 0, 0}, XA4[4] = {0, 0, 0, 0};
//float m1,m2;

String inputString = "";         // a string to hold incoming data
boolean stringComplete = false;  // whether the string is complete
float mov, u1, u2, u3, u4, u5, u6, y2, y3, MAV1, MAV2;
float f1, f2, f3, f4, f5;
long tic = 0, tic2 = 0, tic3 = 0, toc = 0, tec = 0;

// parâmetros de controle pid/GMV
float y = 0, yr = 0, e0 = 0, e1 = 0, e2 = 0, Kp = 0, Kd = 0, Ki = 0, s[3] = { 15.948911710336251 ,-23.381104567352942   ,7.435382639358760};
float U1 = 0, U0 = 0, lei = 0, r4;

const float F[3] = { 12.350837228484437 ,-18.106327376958163   ,5.757960315919436}, lam = 4;
const float BECQ[7] = {1.101132049155242   ,0.064637469732694   ,0.266953082987383   ,0.255304156350001   ,0.176987497484182, 0, 0};
//const float F[3] = {14.740586189276470 , -21.612647470717199   , 6.872061281440709}, lam = 5;
//const float BECQ[7] = {1.120699970333423, 0.123716857428111, 0.315683839632232, 0.304702657118536, 0.304702657118536, 0.304702657118535, 0.211232600092332};
float H[5] = {  15.948911710336251 ,-27.220007616030877  ,13.304043075546691  ,-2.142751280260683   ,0.112274277854317};
float DU[7] = {0, 0, 0, 0, 0, 0, 0},YR[5], emg[5] = {0, 0, 0, 0, 0}, yM[5] = {0, 0, 0, 0, 0}, Z[4];///////////////////////////////////////////////////////////
int mode = 0, refs;
void setup() {
  /////////////////////fast ADC DUE
  ADC->ADC_MR |= 0x80; // these lines set free running mode on adc 7 and adc 6 (pin A0 and A1 - see Due Pinout Diagram thread)
  ADC->ADC_CR = 2;
  ADC->ADC_CHER = 0xBC; // this is (1<<7) | (1<<6) for adc 7 and adc 6
  //
  pinMode(10, OUTPUT);
  pinMode(11, OUTPUT);
  //////////////////////////////////////////////////////////////////
  Serial.begin(250000);
  int t = analogRead(0);

  //Serial.begin(9600);
  // Serial.println("EMG KALMAN filter and xa4zzy classification");
  //ADCSRA &= ~PS_128; //removeo pré_escalamento de 128
  //ADCSRA |= PS_16; // Acrescenta o pé_escalamento de 16 (1MH)

  FK2[0][0] = 1.245124908730679;  FK2[0][1] = 1.707794133368496 ; FK2[0][2] =  1.000000000000000; FK2[0][3] = 0;
  FK2[1][0] = 0.003429842653785;  FK2[1][1] = -0.536234644236777 ; FK2[1][2] = 0; FK2[1][3] =  1.000000000000000;
  FK2[2][0] = 0.003029621643510;  FK2[2][1] = -0.056504970375728; FK2[2][2] =  0; FK2[2][3] = 0;
  FK2[3][0] = -0.300204851316695; FK2[3][1] = -0.115054518755990 ; FK2[3][2] = 0; FK2[3][3] =  0;
  
  FK[0][0] = 1.245124908730679;  FK[0][1] = 1.707794133368496 ; FK[0][2] =  1.000000000000000; FK[0][3] = 0;
  FK[1][0] = 0.003429842653785;  FK[1][1] = -0.536234644236777 ; FK[1][2] = 0; FK[1][3] =  1.000000000000000;
  FK[2][0] = 0.003029621643510;  FK[2][1] = -0.056504970375728; FK[2][2] =  0; FK[2][3] = 0;
  FK[3][0] = -0.300204851316695; FK[3][1] = -0.115054518755990 ; FK[3][2] = 0; FK[3][3] =  0;

  //////////////fusão
  P_[0][0] =  998.5363489356368 ; P_[0][1] =  -190.9128425749228; P_[0][2] = 013.4176443170845; P_[0][3] = -286.2870529166633;
  P_[1][0] = -190.9128425749228 ; P_[1][1] =  044.1877566960686; P_[1][2] = -002.5620016239773 ; P_[1][3] = 055.6017079444693;
  P_[2][0] =  013.4176443170845;  P_[2][1] =  -002.5620016239774; P_[2][2] = 002.2346916949478; P_[2][3] = -004.1765422029379;
  P_[3][0] = -286.2870529166633;  P_[3][1] =  055.6017079444693 ; P_[3][2] = -004.1765422029379 ; P_[3][3] = 087.1050390253165;

  P_2[0][0] =  998.5363489356368 ; P_2[0][1] =  -190.9128425749228; P_2[0][2] = 013.4176443170845; P_2[0][3] = -286.2870529166633;
  P_2[1][0] = -190.9128425749228 ; P_2[1][1] =  044.1877566960686; P_2[1][2] = -002.5620016239773 ; P_2[1][3] = 055.6017079444693;
  P_2[2][0] =  013.4176443170845;  P_2[2][1] =  -002.5620016239774; P_2[2][2] = 002.2346916949478; P_2[2][3] = -004.1765422029379;
  P_2[3][0] = -286.2870529166633;  P_2[3][1] =  055.6017079444693 ; P_2[3][2] = -004.1765422029379 ; P_2[3][3] = 087.1050390253165;

  Copy((float*)P_, 4, 4, (float*)Pa);// Pa recebe P
  Copy((float*)P_, 4, 4, (float*)Pi);// Pi recebe P

  Copy((float*)P_2, 4, 4, (float*)Pa1);// Pa1 recebe P1
  Copy((float*)P_2, 4, 4, (float*)Pi1);// Pi1 recebe P1
  //Print((float*)Pa, 4, 4, "Pa");

  Multiply2((float*)P_, (float)2, 4, 4, 4, (float*)P_1); //multiplica P_1 por 2

  Multiply2((float*)P_, (float)2, 4, 4, 4, (float*)P_); //multiplica P_ por 2

  //Print((float*)P_, 4, 4, "P_");

  Copy((float*)P_, 4, 4, (float*)invP_);// invP recebe P
  Invert((float*)invP_, 4); // inverte e adciona em invP


  Copy((float*)P_, 4, 4, (float*)invP_1);// invP1 recebe P1
  Invert((float*)invP_1, 4); // inverte e adciona em invP

//  Print((float*)invP_, 4, 4, "invP_");

  //////////////////
  for (int z = 0 ; z < 500; z++) {
    gg1 = "";
    analogread();

    f1 = filtroK1(a1); //Aquisição do sinal
    f2 = filtroK2(a2); //Aquisição do sinal
    f3 = filtroK3(a3); //Aquisição do sinal
    f4 = filtroK4(a4); //Aquisição do sinal

    u5 = fuK1k2();
    u6 = fuK3k4();

    mov = u5 - u6;
    //max2=max2+mov;
    max1 = max(max1, mov);
    //max2 = max(u6, max2);
    mix1 = min(mix1, mov);
    //max2=max2+mov;
    max2 = max(a1, max2);
    max3 = max(a2, max3);
    max4 = max(a3, max4);
    max5 = max(a4, max5);
    mix2 = min(a1, mix2);
    mix3 = min(a2, mix3);
    mix4 = min(a3, mix4);
    mix5 = min(a4, mix5);
    mix6 = (mix5, mix4, mix3, mix2);
    max6 = (max5, max4, max3, max2);

    mov = mapfloat(mov, mix1, max1, 0.0, 5.00);
    gg1 = gg1 + max1 + ',' + mix1 + ',' + mov + ',' + u5 + ',' + u6 + ',' + 0.00 + ','  + mov + ',' + 0.00 + ',' + 0.00 + ','  + 0.00 + ',' + 0.00 + ',' + 0.00 + ",E\r\n";
    Serial.print(gg1);
    delayMicroseconds(5000);
  }
  //max2 = mov / 500;
  // Serial.println(max1);
  // Serial.println(mix1);

  toc = micros();
  tic = micros();
  tic2 = micros();
  U0 = 0;  e0 = 0;  i = 0;

  Serial.setTimeout(1);
}

void loop() {

  tec = micros();

  if ((micros() - tic) > 20000) { // tempo de amostragem
    a0 = 0; a1 = 0; a2 = 0; a3 = 0, a4 = 0, i = 0;
    analogread();//fast analog read
    gg1 = ""; gg2 = "";

    u1 = mapfloat(a1, mix2, max2, 0, max6); // conversão para voltagem
    u2 = mapfloat(a2, mix3, max3, 0, max6); // conversão para voltagem
    u3 = mapfloat(a3, mix4, max4, 0, max6); // conversão para voltagem
    u4 = mapfloat(a4, mix5, max5, 0, max6); // conversão para voltagem

    f1 = filtroK1(a1); //Aquisição do sinal
    f2 = filtroK2(a2); //Aquisição do sinal
    f3 = filtroK3(u3); //Aquisição do sinal
    f4 = filtroK4(u4); //Aquisição do sinal

    u5 = fuK1k2();
    u6 = fuK3k4();

    mov = u5 - u6;
    mov = mapfloat(mov, mix1, max1, 0.0, 5.00);
    y = 5.00 * a0 / 4096;

   gg1 = gg1 + a1 + ',' + f1 + ',' + a2 + ',' + f2 + ',' + u5 + ',' + u3 + ','  + f3 + ',' + u4 + ',' + f4 + ','  + u6 + ',' + mov + ',' + r4 + ",E\r\n";

    if (mode == 3) {
      //   gg1 = gg1 + f1 + ',' + f2 + ',' + u5 + ',' + f3 + ',' + f4 + ',' + u6 + ",T\r\n";
    }
    if (mode == 1 || mode == 2) {
      // gg1 = gg1 + U0 + ',' + y + ',' + yr + ',' + e0 + ",T\r\n";
    }
    r4 = micros() - tec;
    //gg1 = gg1 + U0 + ',' + y + ',' + yr + ',' + e0 + ",T\r\n";

    Serial.print(gg1);
    tic = micros();
  }

  if (mode == 1) {
    refRand(5);
    //yr = mov;
    GMV();
  }
  if (mode == 2) {
    refRand(5);
    //yr = mov;
    PID();
  }

  if (Serial.available() > 0) {
   
    char READ = Serial.read();
    if (READ == 'p') {
      mode = 0;
      Serial.print("Parado");
    }
    if (READ == 's') {
      mode = 1;
      Serial.print("Gmv");
    }
    if (READ == 'd') {
      mode = 2;
      Serial.print("PID");
    }
    U0 = Serial.parseFloat();
    saida(U0);
    u1 = 5.00 * a0 / 4906;
    Serial.println(u1);

  }
}

///////////////////////funções//////////////////////////////////////
void refRand(long T) {
  T = T * 1000000;

  if ((micros() - tic3) > T) { // tempo de amostragem
    if (refs == 1) {
      yr = 4;
      refs = 0;
    }
    else if (refs == 0) {
      yr = 1;
      refs = 1;
    }
    tic3 = micros();
  }

}
void GMV() {
  if ((micros() - tic2) > 100000) { // tempo de amostragem
    emg[0] = yr; // recebe o valor de referencia angular
    yM[0] = y; //recebe valor de posição lido do sensor do motor

    Z[0] = (-BECQ[1] * DU[1] - BECQ[2] * DU[2] - BECQ[3] * DU[3] - BECQ[4] * DU[4] - BECQ[5] * DU[5] - BECQ[6] * DU[6]);

    Z[1] = (H[0] * emg[0] + H[1] * emg[1] + H[2] * emg[2] + H[3] * emg[3] + H[4] * emg[4]);

    Z[2] = (-F[0] * yM[0] - F[1] * yM[1] - F[2] * yM[2]);

    Z[3] = Z[0] + Z[1] + Z[2];
    DU[0] = (Z[3] / BECQ[0]);

    U0 = U1 + DU[0];
    e0 = yr - y;

    saida(U0);

    if (i > 500000) {
      // U0 = 0;    e0 = 0;    i = 0;
    }
    U1 = U0; DU[6] = DU[5]; DU[5] = DU[4]; DU[4] = DU[3]; DU[3] = DU[2]; DU[2] = DU[1]; DU[1] = DU[0]; emg[4] = emg[3]; emg[3] = emg[2]; emg[2] = emg[1]; emg[1] = emg[0]; yM[2] = yM[1]; yM[1] = yM[0];
    i++;

    gg1 = gg1 + U0 + ',' + y + ',' +yr + ',' + e0 + ",T\r\n";
    Serial.print(gg1);

    tic2 = micros();
  }
}
void PID() {
  if ((micros() - tic2) > 100000) { // tempo de amostragem
    // yr=100;
    //y=(analogRead(0));
    e0 = yr - y;
    U0 = U1 + s[0] * e0 + s[1] * e1 + s[2] * e2;
    saida(U0);

    if (i > 500) {
      // U0 = 0;    e0 = 0;    i = 0;
    }
    e2 = e1;   e1 = e0;    U1 = U0;
    i++;

    gg1 = gg1 + U0 + ',' + y + ',' + yr + ',' + e0 + ",T\r\n";
    Serial.print(gg1);

    tic2 = micros();
  }
}
void saida(float lei) {

  if (lei < 0) {
    lei = mapfloat(-lei, 0, 5, 85, 255);
    lei=-lei;
  }
  if (lei > 0) {
    lei = mapfloat(lei, 0, 5, 85, 255);
  }
  
  if (lei < -255) {
    lei = -255;
  }
  if (lei > 255) {
    lei = 255;
  }
  // Serial.print(u[0]);
  if (lei < 0) {
    analogWrite(11, -lei);
    analogWrite(10, 0);
  }
  else if (lei >= 0) {
    analogWrite(11, 0);
    analogWrite(10, lei);
    //  Serial.print("lei1: ");
    // Serial.println(lei);
  }

}
float mapfloat(long x, long in_min, long in_max, long out_min, long out_max)
{
  return (float)(x - in_min) * (out_max - out_min) / (float)(in_max - in_min) + out_min;
}


int analogread() {
while ((ADC->ADC_ISR & 0xBC) != 0xBC); // wait for two conversions (pin A0[7]  and A1[6])
  a0 = ADC->ADC_CDR[7];            // read data on A0 pin
  a1 = ADC->ADC_CDR[5];            // read data on A0 pin
  a2 = ADC->ADC_CDR[4];            // read data on A1 pin
  a3 = ADC->ADC_CDR[3];            // read data on A0 pin
  a4 = ADC->ADC_CDR[2];            // read data on A1 pin


}


float filtroK1(float u) {
  //Filtro de kalman
  Multiply((float*)FK, (float*)xa1, 4, 4, 1, (float*)C);//C=(Ad-L*Cd)*xa
  for (int k = 0; k < 4; k++) {
    L1[k] = L_1[k] * u;
  }
  Add((float*) C, (float*) L1, 4, 1, (float*) xa1);//xa=C+l
  Multiply((float*)Cd, (float*)xa1, 4, 4, 1, (float*)B);//B=Cd*xa
  Copy((float*)xa1, 4, 1, (float*)XA1);
  return B[0];
}
float filtroK2(float u ) {
  //Filtro de kalman
  Multiply((float*)FK, (float*)xa2, 4, 4, 1, (float*)C);//C=(Ad-L*Cd)*xa
  for (int k = 0; k < 4; k++) {
    L2[k] = L_1[k] * u;
  }
  Add((float*) C, (float*) L2, 4, 1, (float*) xa2);//xa=C+l
  Multiply((float*)Cd, (float*)xa2, 4, 4, 1, (float*)B2);//B=Cd*xa
  Copy((float*)xa2, 4, 1, (float*)XA2);

  return B2[0];
}
float filtroK3(float u ) {
  //Filtro de kalman
  Multiply((float*)FK2, (float*)xa3, 4, 4, 1, (float*)C);//C=(Ad-L*Cd)*xa
  for (int k = 0; k < 4; k++) {
    L3[k] = L_2[k] * u;
  }
  Add((float*) C, (float*) L3, 4, 1, (float*) xa3);//xa=C+l
  Multiply((float*)Cd, (float*)xa3, 4, 4, 1, (float*)B3);//B=Cd*xa
  Copy((float*)xa3, 4, 1, (float*)XA3);

  return B3[0];
}
float filtroK4(float u ) {
  //Filtro de kalman
  Multiply((float*)FK2, (float*)xa4, 4, 4, 1, (float*)C);//C=(Ad-L*Cd)*xa
  for (int k = 0; k < 4; k++) {
    L4[k] = L_2[k] * u;
  }
  Add((float*) C, (float*) L4, 4, 1, (float*) xa4);//xa=C+l
  Multiply((float*)Cd, (float*)xa4, 4, 4, 1, (float*)B4);//B=Cd*xa
  Copy((float*)xa4, 4, 1, (float*)XA4);
  return B4[0];
}


void Add(float* A, float* B, int m, int n, float* C)
{
  int i, j;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C[n * i + j] = A[n * i + j] + B[n * i + j];
}
void Multiply(float* A, float* B, int m, int p, int n, float* C)
{
  int i, j, k;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      C[n * i + j] = 0;
      for (k = 0; k < p; k++)
        C[n * i + j] = C[n * i + j] + A[p * i + k] * B[n * k + j];
    }
}
void sub(float* A, float* B, int m, int n, float* C)
{
  int i, j;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      C[n * i + j] = A[n * i + j] - B[n * i + j];
}

void Multiply2(float* A, float B, int m, int p, int n, float* C)
{
  int i, j, k;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
    {
      C[n * i + j] = A[p * i + j] * B;
    }
  }
}
void Copy(float* A, int n, int m, float* B)
{
  int i, j, k;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      B[n * i + j] = A[n * i + j];
    }
}

float fuK1k2() {

  //x1= inv(P1_)*( P1_*x1 + Pa*(xa1-m1)+Pi*(xa2-m2) );
  sub((float*) XA1, (float*) m1, 4, 1, (float*) M1);//(xa1-m1)
  sub((float*) XA2, (float*) m2, 4, 1, (float*) M2);//(xa2-m2)
  //Print((float*)M2, 4, 1, "m2");
  //Print((float*)xa1, 4, 1, "xa1");
  Multiply((float*)Pa, (float*)M1, 4, 4, 1, (float*)m1);//Pa*(xa1-m1)
  Multiply((float*)Pi, (float*)M2, 4, 4, 1, (float*)m2);//Pi*(xa2-m2)

  //Print((float*)Pa, 4, 4, "Pa");
  //Print((float*)m1, 4, 1, "m1");

  Add((float*) m1, (float*) m2, 4, 1, (float*) M1);//Pa*(xa1-m1)+Pi*(xa2-m2)

  // Print((float*)M1, 4, 1, "m1");
  //Print((float*)xa4, 4, 1, "xa4");

  Multiply((float*)P_, (float*)xa12, 4, 4, 1, (float*)M2);//P1_*x1

  //Print((float*)M2, 4, 1, "M2");

  Add((float*) M2, (float*) M1, 4, 1, (float*) m2);//( P1_*x1 + Pa*(xa1-m1)+Pi*(xa2-m2) )

  Multiply((float*)invP_, (float*)m2, 4, 4, 1, (float*)xa12);//Pa*(xa1-m1)

  //Print((float*)xa4, 4, 1, "xa4");
  // y1(i)=Cd*x1;
  Multiply((float*)Cd, (float*)xa12, 4, 4, 1, (float*)B12);//B=Cd*xa
  //Print((float*)B4, 4, 1, "B4");

  Copy((float*)XA1, 4, 1, (float*)m1);
  Copy((float*)XA2, 4, 1, (float*)m2);
  //Serial.println(B4[0]);
  return B12[0];
}

float fuK3k4() {

  //x1= inv(P1_)*( P1_*x1 + Pa*(xa1-m1)+Pi*(xa2-m2) );
  sub((float*) XA3, (float*) m3, 4, 1, (float*) M3);//(xa1-m1)
  sub((float*) XA4, (float*) m4, 4, 1, (float*) M4);//(xa2-m2)
  //Print((float*)M2, 4, 1, "m2");
  //Print((float*)xa1, 4, 1, "xa1");
  Multiply((float*)Pa1, (float*)M3, 4, 4, 1, (float*)m3);//Pa*(xa1-m1)
  Multiply((float*)Pi1, (float*)M4, 4, 4, 1, (float*)m4);//Pi*(xa2-m2)

  //Print((float*)Pa, 4, 4, "Pa");
  //Print((float*)m1, 4, 1, "m1");

  Add((float*) m3, (float*) m4, 4, 1, (float*) M3);//Pa*(xa1-m1)+Pi*(xa2-m2)

  // Print((float*)M1, 4, 1, "m1");
  //Print((float*)xa4, 4, 1, "xa4");

  Multiply((float*)P_1, (float*)xa34, 4, 4, 1, (float*)M4);//P1_*x1

  //Print((float*)M2, 4, 1, "M2");

  Add((float*) M4, (float*) M3, 4, 1, (float*) m4);//( P1_*x1 + Pa*(xa1-m1)+Pi*(xa2-m2) )

  Multiply((float*)invP_1, (float*)m4, 4, 4, 1, (float*)xa34);//Pa*(xa1-m1)

  //Print((float*)xa4, 4, 1, "xa4");
  // y1(i)=Cd*x1;
  Multiply((float*)Cd, (float*)xa34, 4, 4, 1, (float*)B34);//B=Cd*xa
  //Print((float*)B4, 4, 1, "B4");

  Copy((float*)XA3, 4, 1, (float*)m3);
  Copy((float*)XA4, 4, 1, (float*)m4);
  //Serial.println(B4[0]);
  return B34[0];
}
int Invert(float* A, int n)
{
  // A = input matrix AND result matrix
  // n = number of rows = number of columns in A (n x n)
  int pivrow;   // keeps track of current pivot row
  int k, i, j;    // k: overall index along diagonal; i: row index; j: col index
  int pivrows[n]; // keeps track of rows swaps to undo at end
  float tmp;    // used for finding max value and making column swaps

  for (k = 0; k < n; k++)
  {
    // find pivot row, the row with biggest entry in current column
    tmp = 0;
    for (i = k; i < n; i++)
    {
      if (abs(A[i * n + k]) >= tmp) // 'Avoid using other functions inside abs()?'
      {
        tmp = abs(A[i * n + k]);
        pivrow = i;
      }
    }

    // check for singular matrix
    if (A[pivrow * n + k] == 0.0f)
    {
      Serial.println("Inversion failed due to singular matrix");
      return 0;
    }

    // Execute pivot (row swap) if needed
    if (pivrow != k)
    {
      // swap row k with pivrow
      for (j = 0; j < n; j++)
      {
        tmp = A[k * n + j];
        A[k * n + j] = A[pivrow * n + j];
        A[pivrow * n + j] = tmp;
      }
    }
    pivrows[k] = pivrow;  // record row swap (even if no swap happened)

    tmp = 1.0f / A[k * n + k];  // invert pivot element
    A[k * n + k] = 1.0f;    // This element of input matrix becomes result matrix

    // Perform row reduction (divide every element by pivot)
    for (j = 0; j < n; j++)
    {
      A[k * n + j] = A[k * n + j] * tmp;
    }

    // Now eliminate all other entries in this column
    for (i = 0; i < n; i++)
    {
      if (i != k)
      {
        tmp = A[i * n + k];
        A[i * n + k] = 0.0f; // The other place where in matrix becomes result mat
        for (j = 0; j < n; j++)
        {
          A[i * n + j] = A[i * n + j] - A[k * n + j] * tmp;
        }
      }
    }
  }

  // Done, now need to undo pivot row swaps by doing column swaps in reverse order
  for (k = n - 1; k >= 0; k--)
  {
    if (pivrows[k] != k)
    {
      for (i = 0; i < n; i++)
      {
        tmp = A[i * n + k];
        A[i * n + k] = A[i * n + pivrows[k]];
        A[i * n + pivrows[k]] = tmp;
      }
    }
  }
  return 1;
}

void normali() {

}
