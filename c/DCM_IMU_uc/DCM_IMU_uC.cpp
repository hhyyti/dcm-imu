#include "DCM_IMU_uC.h"

DCM_IMU_uC::DCM_IMU_uC(const float Gravity, const float *State, const float *Covariance,
		const float DCMVariance, const float BiasVariance,
		const float InitialDCMVariance, const float InitialBiasVariance,
		const float MeasurementVariance, const float MeasurementVarianceVariableGain) :
		g0(Gravity), q_dcm2(DCMVariance), q_gyro_bias2(BiasVariance),
		r_acc2(MeasurementVariance), r_a2(MeasurementVarianceVariableGain) {

	float temp[] = DEFAULT_state;
	if (State != NULL) {
		for (int i = 0; i < 6; ++i) {
			temp[i] = State[i];
		}
	}
	x0 = temp[0];
	x1 = temp[1];
	x2 = temp[2];
	x3 = temp[3];
	x4 = temp[4];
	x5 = temp[5];

	if (Covariance == NULL) {
		P00 = InitialDCMVariance;
		P01 = 0;
		P02 = 0;
		P03 = 0;
		P04 = 0;
		P05 = 0;
		P10 = 0;
		P11 = InitialDCMVariance;
		P12 = 0;
		P13 = 0;
		P14 = 0;
		P15 = 0;
		P20 = 0;
		P21 = 0;
		P22 = InitialDCMVariance;
		P23 = 0;
		P24 = 0;
		P25 = 0;
		P30 = 0;
		P31 = 0;
		P32 = 0;
		P33 = InitialBiasVariance;
		P34 = 0;
		P35 = 0;
		P40 = 0;
		P41 = 0;
		P42 = 0;
		P43 = 0;
		P44 = InitialBiasVariance;
		P45 = 0;
		P50 = 0;
		P51 = 0;
		P52 = 0;
		P53 = 0;
		P54 = 0;
		P55 = InitialBiasVariance;
	}
	else {
		P00 = Covariance[0];
		P01 = Covariance[1];
		P02 = Covariance[2];
		P03 = Covariance[3];
		P04 = Covariance[4];
		P05 = Covariance[5];
		P10 = Covariance[6];
		P11 = Covariance[7];
		P12 = Covariance[8];
		P13 = Covariance[9];
		P14 = Covariance[10];
		P15 = Covariance[11];
		P20 = Covariance[12];
		P21 = Covariance[13];
		P22 = Covariance[14];
		P23 = Covariance[15];
		P24 = Covariance[16];
		P25 = Covariance[17];
		P30 = Covariance[18];
		P31 = Covariance[19];
		P32 = Covariance[20];
		P33 = Covariance[21];
		P34 = Covariance[22];
		P35 = Covariance[23];
		P40 = Covariance[24];
		P41 = Covariance[25];
		P42 = Covariance[26];
		P43 = Covariance[27];
		P44 = Covariance[28];
		P45 = Covariance[29];
		P50 = Covariance[30];
		P51 = Covariance[31];
		P52 = Covariance[32];
		P53 = Covariance[33];
		P54 = Covariance[34];
		P55 = Covariance[35];
	}

	g0_2 = g0*g0;
	yaw = 0.0f;
	pitch = 0.0f;
	roll = 0.0f;
}


void DCM_IMU_uC::updateIMU(const float *Gyroscope, const float *Accelerometer, const float h) {

	// save last state to memory for rotation estimation
	float x_last[3];
	x_last[0] = x0;
	x_last[1] = x1;
	x_last[2] = x2;

	// control input (gyroscopes)
	float u0 = Gyroscope[0];
	float u1 = Gyroscope[1];
	float u2 = Gyroscope[2];

	// state prediction
	float x_0 = x0-h*(u1*x2-u2*x1+x1*x5-x2*x4);
	float x_1 = x1+h*(u0*x2-u2*x0+x0*x5-x2*x3);
	float x_2 = x2-h*(u0*x1-u1*x0+x0*x4-x1*x3);
	float x_3 = x3;
	float x_4 = x4;
	float x_5 = x5;

	// covariance prediction
	float hh = h*h;
	float P_00 = P00-h*(P05*x1-P04*x2-P40*x2+P50*x1+P02*(u1-x4)+P20*(u1-x4)-P01*(u2-x5)-P10*(u2-x5))+hh*(q_dcm2-x1*(P45*x2-P55*x1-P25*(u1-x4)+P15*(u2-x5))+x2*(P44*x2-P54*x1-P24*(u1-x4)+P14*(u2-x5))-(u1-x4)*(P42*x2-P52*x1-P22*(u1-x4)+P12*(u2-x5))+(u2-x5)*(P41*x2-P51*x1-P21*(u1-x4)+P11*(u2-x5)));
	float P_01 = P01+h*(P05*x0-P03*x2+P41*x2-P51*x1+P02*(u0-x3)-P00*(u2-x5)-P21*(u1-x4)+P11*(u2-x5))+hh*(x0*(P45*x2-P55*x1-P25*(u1-x4)+P15*(u2-x5))-x2*(P43*x2-P53*x1-P23*(u1-x4)+P13*(u2-x5))+(u0-x3)*(P42*x2-P52*x1-P22*(u1-x4)+P12*(u2-x5))-(u2-x5)*(P40*x2-P50*x1-P20*(u1-x4)+P10*(u2-x5)));
	float P_02 = P02-h*(P04*x0-P03*x1-P42*x2+P52*x1+P01*(u0-x3)-P00*(u1-x4)+P22*(u1-x4)-P12*(u2-x5))-hh*(x0*(P44*x2-P54*x1-P24*(u1-x4)+P14*(u2-x5))-x1*(P43*x2-P53*x1-P23*(u1-x4)+P13*(u2-x5))+(u0-x3)*(P41*x2-P51*x1-P21*(u1-x4)+P11*(u2-x5))-(u1-x4)*(P40*x2-P50*x1-P20*(u1-x4)+P10*(u2-x5)));
	float P_03 = P03+h*(P43*x2-P53*x1-P23*(u1-x4)+P13*(u2-x5));
	float P_04 = P04+h*(P44*x2-P54*x1-P24*(u1-x4)+P14*(u2-x5));
	float P_05 = P05+h*(P45*x2-P55*x1-P25*(u1-x4)+P15*(u2-x5));
	float P_10 = P10-h*(P15*x1-P14*x2+P30*x2-P50*x0-P20*(u0-x3)+P12*(u1-x4)+P00*(u2-x5)-P11*(u2-x5))+hh*(x1*(P35*x2-P55*x0-P25*(u0-x3)+P05*(u2-x5))-x2*(P34*x2-P54*x0-P24*(u0-x3)+P04*(u2-x5))+(u1-x4)*(P32*x2-P52*x0-P22*(u0-x3)+P02*(u2-x5))-(u2-x5)*(P31*x2-P51*x0-P21*(u0-x3)+P01*(u2-x5)));
	float P_11 = P11+h*(P15*x0-P13*x2-P31*x2+P51*x0+P12*(u0-x3)+P21*(u0-x3)-P01*(u2-x5)-P10*(u2-x5))+hh*(q_dcm2-x0*(P35*x2-P55*x0-P25*(u0-x3)+P05*(u2-x5))+x2*(P33*x2-P53*x0-P23*(u0-x3)+P03*(u2-x5))-(u0-x3)*(P32*x2-P52*x0-P22*(u0-x3)+P02*(u2-x5))+(u2-x5)*(P30*x2-P50*x0-P20*(u0-x3)+P00*(u2-x5)));
	float P_12 = P12-h*(P14*x0-P13*x1+P32*x2-P52*x0+P11*(u0-x3)-P22*(u0-x3)-P10*(u1-x4)+P02*(u2-x5))+hh*(x0*(P34*x2-P54*x0-P24*(u0-x3)+P04*(u2-x5))-x1*(P33*x2-P53*x0-P23*(u0-x3)+P03*(u2-x5))+(u0-x3)*(P31*x2-P51*x0-P21*(u0-x3)+P01*(u2-x5))-(u1-x4)*(P30*x2-P50*x0-P20*(u0-x3)+P00*(u2-x5)));
	float P_13 = P13-h*(P33*x2-P53*x0-P23*(u0-x3)+P03*(u2-x5));
	float P_14 = P14-h*(P34*x2-P54*x0-P24*(u0-x3)+P04*(u2-x5));
	float P_15 = P15-h*(P35*x2-P55*x0-P25*(u0-x3)+P05*(u2-x5));
	float P_20 = P20-h*(P25*x1-P30*x1+P40*x0-P24*x2+P10*(u0-x3)-P00*(u1-x4)+P22*(u1-x4)-P21*(u2-x5))-hh*(x1*(P35*x1-P45*x0-P15*(u0-x3)+P05*(u1-x4))-x2*(P34*x1-P44*x0-P14*(u0-x3)+P04*(u1-x4))+(u1-x4)*(P32*x1-P42*x0-P12*(u0-x3)+P02*(u1-x4))-(u2-x5)*(P31*x1-P41*x0-P11*(u0-x3)+P01*(u1-x4)));
	float P_21 = P21+h*(P25*x0+P31*x1-P41*x0-P23*x2-P11*(u0-x3)+P01*(u1-x4)+P22*(u0-x3)-P20*(u2-x5))+hh*(x0*(P35*x1-P45*x0-P15*(u0-x3)+P05*(u1-x4))-x2*(P33*x1-P43*x0-P13*(u0-x3)+P03*(u1-x4))+(u0-x3)*(P32*x1-P42*x0-P12*(u0-x3)+P02*(u1-x4))-(u2-x5)*(P30*x1-P40*x0-P10*(u0-x3)+P00*(u1-x4)));
	float P_22 = P22-h*(P24*x0-P23*x1-P32*x1+P42*x0+P12*(u0-x3)+P21*(u0-x3)-P02*(u1-x4)-P20*(u1-x4))+hh*(q_dcm2-x0*(P34*x1-P44*x0-P14*(u0-x3)+P04*(u1-x4))+x1*(P33*x1-P43*x0-P13*(u0-x3)+P03*(u1-x4))-(u0-x3)*(P31*x1-P41*x0-P11*(u0-x3)+P01*(u1-x4))+(u1-x4)*(P30*x1-P40*x0-P10*(u0-x3)+P00*(u1-x4)));
	float P_23 = P23+h*(P33*x1-P43*x0-P13*(u0-x3)+P03*(u1-x4));
	float P_24 = P24+h*(P34*x1-P44*x0-P14*(u0-x3)+P04*(u1-x4));
	float P_25 = P25+h*(P35*x1-P45*x0-P15*(u0-x3)+P05*(u1-x4));
	float P_30 = P30-h*(P35*x1-P34*x2+P32*(u1-x4)-P31*(u2-x5));
	float P_31 = P31+h*(P35*x0-P33*x2+P32*(u0-x3)-P30*(u2-x5));
	float P_32 = P32-h*(P34*x0-P33*x1+P31*(u0-x3)-P30*(u1-x4));
	float P_33 = P33+hh*q_gyro_bias2;
	float P_34 = P34;
	float P_35 = P35;
	float P_40 = P40-h*(P45*x1-P44*x2+P42*(u1-x4)-P41*(u2-x5));
	float P_41 = P41+h*(P45*x0-P43*x2+P42*(u0-x3)-P40*(u2-x5));
	float P_42 = P42-h*(P44*x0-P43*x1+P41*(u0-x3)-P40*(u1-x4));
	float P_43 = P43;
	float P_44 = P44+hh*q_gyro_bias2;
	float P_45 = P45;
	float P_50 = P50-h*(P55*x1-P54*x2+P52*(u1-x4)-P51*(u2-x5));
	float P_51 = P51+h*(P55*x0-P53*x2+P52*(u0-x3)-P50*(u2-x5));
	float P_52 = P52-h*(P54*x0-P53*x1+P51*(u0-x3)-P50*(u1-x4));
	float P_53 = P53;
	float P_54 = P54;
	float P_55 = P55+hh*q_gyro_bias2;

	// measurements (accelerometers)
	float z0 = Accelerometer[0];
	float z1 = Accelerometer[1];
	float z2 = Accelerometer[2];

	// Kalman innovation
	float y0 = z0-g0*x_0;
	float y1 = z1-g0*x_1;
	float y2 = z2-g0*x_2;

	float a_len = sqrt(y0*y0+y1*y1+y2*y2);

	float S00 = r_acc2+a_len*r_a2+P_00*g0_2;
	float S01 = P_01*g0_2;
	float S02 = P_02*g0_2;
	float S10 = P_10*g0_2;
	float S11 = r_acc2+a_len*r_a2+P_11*g0_2;
	float S12 = P_12*g0_2;
	float S20 = P_20*g0_2;
	float S21 = P_21*g0_2;
	float S22 = r_acc2+a_len*r_a2+P_22*g0_2;

	// Kalman gain
	float invPart = 1.0 / (S00*S11*S22-S00*S12*S21-S01*S10*S22+S01*S12*S20+S02*S10*S21-S02*S11*S20);
	float K00 = (g0*(P_02*S10*S21-P_02*S11*S20-P_01*S10*S22+P_01*S12*S20+P_00*S11*S22-P_00*S12*S21))*invPart;
	float K01 = -(g0*(P_02*S00*S21-P_02*S01*S20-P_01*S00*S22+P_01*S02*S20+P_00*S01*S22-P_00*S02*S21))*invPart;
	float K02 = (g0*(P_02*S00*S11-P_02*S01*S10-P_01*S00*S12+P_01*S02*S10+P_00*S01*S12-P_00*S02*S11))*invPart;
	float K10 = (g0*(P_12*S10*S21-P_12*S11*S20-P_11*S10*S22+P_11*S12*S20+P_10*S11*S22-P_10*S12*S21))*invPart;
	float K11 = -(g0*(P_12*S00*S21-P_12*S01*S20-P_11*S00*S22+P_11*S02*S20+P_10*S01*S22-P_10*S02*S21))*invPart;
	float K12 = (g0*(P_12*S00*S11-P_12*S01*S10-P_11*S00*S12+P_11*S02*S10+P_10*S01*S12-P_10*S02*S11))*invPart;
	float K20 = (g0*(P_22*S10*S21-P_22*S11*S20-P_21*S10*S22+P_21*S12*S20+P_20*S11*S22-P_20*S12*S21))*invPart;
	float K21 = -(g0*(P_22*S00*S21-P_22*S01*S20-P_21*S00*S22+P_21*S02*S20+P_20*S01*S22-P_20*S02*S21))*invPart;
	float K22 = (g0*(P_22*S00*S11-P_22*S01*S10-P_21*S00*S12+P_21*S02*S10+P_20*S01*S12-P_20*S02*S11))*invPart;
	float K30 = (g0*(P_32*S10*S21-P_32*S11*S20-P_31*S10*S22+P_31*S12*S20+P_30*S11*S22-P_30*S12*S21))*invPart;
	float K31 = -(g0*(P_32*S00*S21-P_32*S01*S20-P_31*S00*S22+P_31*S02*S20+P_30*S01*S22-P_30*S02*S21))*invPart;
	float K32 = (g0*(P_32*S00*S11-P_32*S01*S10-P_31*S00*S12+P_31*S02*S10+P_30*S01*S12-P_30*S02*S11))*invPart;
	float K40 = (g0*(P_42*S10*S21-P_42*S11*S20-P_41*S10*S22+P_41*S12*S20+P_40*S11*S22-P_40*S12*S21))*invPart;
	float K41 = -(g0*(P_42*S00*S21-P_42*S01*S20-P_41*S00*S22+P_41*S02*S20+P_40*S01*S22-P_40*S02*S21))*invPart;
	float K42 = (g0*(P_42*S00*S11-P_42*S01*S10-P_41*S00*S12+P_41*S02*S10+P_40*S01*S12-P_40*S02*S11))*invPart;
	float K50 = (g0*(P_52*S10*S21-P_52*S11*S20-P_51*S10*S22+P_51*S12*S20+P_50*S11*S22-P_50*S12*S21))*invPart;
	float K51 = -(g0*(P_52*S00*S21-P_52*S01*S20-P_51*S00*S22+P_51*S02*S20+P_50*S01*S22-P_50*S02*S21))*invPart;
	float K52 = (g0*(P_52*S00*S11-P_52*S01*S10-P_51*S00*S12+P_51*S02*S10+P_50*S01*S12-P_50*S02*S11))*invPart;

	// update a posteriori
	x0 = x_0+K00*y0+K01*y1+K02*y2;
	x1 = x_1+K10*y0+K11*y1+K12*y2;
	x2 = x_2+K20*y0+K21*y1+K22*y2;
	x3 = x_3+K30*y0+K31*y1+K32*y2;
	x4 = x_4+K40*y0+K41*y1+K42*y2;
	x5 = x_5+K50*y0+K51*y1+K52*y2;

	// update a posteriori covariance
	float r_adab = (r_acc2+a_len*r_a2);
	float P__00 = P_00-g0*(K00*P_00*2.0+K01*P_01+K01*P_10+K02*P_02+K02*P_20)+(K00*K00)*r_adab+(K01*K01)*r_adab+(K02*K02)*r_adab+g0_2*(K00*(K00*P_00+K01*P_10+K02*P_20)+K01*(K00*P_01+K01*P_11+K02*P_21)+K02*(K00*P_02+K01*P_12+K02*P_22));
	float P__01 = P_01-g0*(K00*P_01+K01*P_11+K02*P_21+K10*P_00+K11*P_01+K12*P_02)+g0_2*(K10*(K00*P_00+K01*P_10+K02*P_20)+K11*(K00*P_01+K01*P_11+K02*P_21)+K12*(K00*P_02+K01*P_12+K02*P_22))+K00*K10*r_adab+K01*K11*r_adab+K02*K12*r_adab;
	float P__02 = P_02-g0*(K00*P_02+K01*P_12+K02*P_22+K20*P_00+K21*P_01+K22*P_02)+g0_2*(K20*(K00*P_00+K01*P_10+K02*P_20)+K21*(K00*P_01+K01*P_11+K02*P_21)+K22*(K00*P_02+K01*P_12+K02*P_22))+K00*K20*r_adab+K01*K21*r_adab+K02*K22*r_adab;
	float P__03 = P_03-g0*(K00*P_03+K01*P_13+K02*P_23+K30*P_00+K31*P_01+K32*P_02)+g0_2*(K30*(K00*P_00+K01*P_10+K02*P_20)+K31*(K00*P_01+K01*P_11+K02*P_21)+K32*(K00*P_02+K01*P_12+K02*P_22))+K00*K30*r_adab+K01*K31*r_adab+K02*K32*r_adab;
	float P__04 = P_04-g0*(K00*P_04+K01*P_14+K02*P_24+K40*P_00+K41*P_01+K42*P_02)+g0_2*(K40*(K00*P_00+K01*P_10+K02*P_20)+K41*(K00*P_01+K01*P_11+K02*P_21)+K42*(K00*P_02+K01*P_12+K02*P_22))+K00*K40*r_adab+K01*K41*r_adab+K02*K42*r_adab;
	float P__05 = P_05-g0*(K00*P_05+K01*P_15+K02*P_25+K50*P_00+K51*P_01+K52*P_02)+g0_2*(K50*(K00*P_00+K01*P_10+K02*P_20)+K51*(K00*P_01+K01*P_11+K02*P_21)+K52*(K00*P_02+K01*P_12+K02*P_22))+K00*K50*r_adab+K01*K51*r_adab+K02*K52*r_adab;
	float P__10 = P_10-g0*(K00*P_10+K01*P_11+K02*P_12+K10*P_00+K11*P_10+K12*P_20)+g0_2*(K00*(K10*P_00+K11*P_10+K12*P_20)+K01*(K10*P_01+K11*P_11+K12*P_21)+K02*(K10*P_02+K11*P_12+K12*P_22))+K00*K10*r_adab+K01*K11*r_adab+K02*K12*r_adab;
	float P__11 = P_11-g0*(K10*P_01+K10*P_10+K11*P_11*2.0+K12*P_12+K12*P_21)+(K10*K10)*r_adab+(K11*K11)*r_adab+(K12*K12)*r_adab+g0_2*(K10*(K10*P_00+K11*P_10+K12*P_20)+K11*(K10*P_01+K11*P_11+K12*P_21)+K12*(K10*P_02+K11*P_12+K12*P_22));
	float P__12 = P_12-g0*(K10*P_02+K11*P_12+K12*P_22+K20*P_10+K21*P_11+K22*P_12)+g0_2*(K20*(K10*P_00+K11*P_10+K12*P_20)+K21*(K10*P_01+K11*P_11+K12*P_21)+K22*(K10*P_02+K11*P_12+K12*P_22))+K10*K20*r_adab+K11*K21*r_adab+K12*K22*r_adab;
	float P__13 = P_13-g0*(K10*P_03+K11*P_13+K12*P_23+K30*P_10+K31*P_11+K32*P_12)+g0_2*(K30*(K10*P_00+K11*P_10+K12*P_20)+K31*(K10*P_01+K11*P_11+K12*P_21)+K32*(K10*P_02+K11*P_12+K12*P_22))+K10*K30*r_adab+K11*K31*r_adab+K12*K32*r_adab;
	float P__14 = P_14-g0*(K10*P_04+K11*P_14+K12*P_24+K40*P_10+K41*P_11+K42*P_12)+g0_2*(K40*(K10*P_00+K11*P_10+K12*P_20)+K41*(K10*P_01+K11*P_11+K12*P_21)+K42*(K10*P_02+K11*P_12+K12*P_22))+K10*K40*r_adab+K11*K41*r_adab+K12*K42*r_adab;
	float P__15 = P_15-g0*(K10*P_05+K11*P_15+K12*P_25+K50*P_10+K51*P_11+K52*P_12)+g0_2*(K50*(K10*P_00+K11*P_10+K12*P_20)+K51*(K10*P_01+K11*P_11+K12*P_21)+K52*(K10*P_02+K11*P_12+K12*P_22))+K10*K50*r_adab+K11*K51*r_adab+K12*K52*r_adab;
	float P__20 = P_20-g0*(K00*P_20+K01*P_21+K02*P_22+K20*P_00+K21*P_10+K22*P_20)+g0_2*(K00*(K20*P_00+K21*P_10+K22*P_20)+K01*(K20*P_01+K21*P_11+K22*P_21)+K02*(K20*P_02+K21*P_12+K22*P_22))+K00*K20*r_adab+K01*K21*r_adab+K02*K22*r_adab;
	float P__21 = P_21-g0*(K10*P_20+K11*P_21+K12*P_22+K20*P_01+K21*P_11+K22*P_21)+g0_2*(K10*(K20*P_00+K21*P_10+K22*P_20)+K11*(K20*P_01+K21*P_11+K22*P_21)+K12*(K20*P_02+K21*P_12+K22*P_22))+K10*K20*r_adab+K11*K21*r_adab+K12*K22*r_adab;
	float P__22 = P_22-g0*(K20*P_02+K20*P_20+K21*P_12+K21*P_21+K22*P_22*2.0)+(K20*K20)*r_adab+(K21*K21)*r_adab+(K22*K22)*r_adab+g0_2*(K20*(K20*P_00+K21*P_10+K22*P_20)+K21*(K20*P_01+K21*P_11+K22*P_21)+K22*(K20*P_02+K21*P_12+K22*P_22));
	float P__23 = P_23-g0*(K20*P_03+K21*P_13+K22*P_23+K30*P_20+K31*P_21+K32*P_22)+g0_2*(K30*(K20*P_00+K21*P_10+K22*P_20)+K31*(K20*P_01+K21*P_11+K22*P_21)+K32*(K20*P_02+K21*P_12+K22*P_22))+K20*K30*r_adab+K21*K31*r_adab+K22*K32*r_adab;
	float P__24 = P_24-g0*(K20*P_04+K21*P_14+K22*P_24+K40*P_20+K41*P_21+K42*P_22)+g0_2*(K40*(K20*P_00+K21*P_10+K22*P_20)+K41*(K20*P_01+K21*P_11+K22*P_21)+K42*(K20*P_02+K21*P_12+K22*P_22))+K20*K40*r_adab+K21*K41*r_adab+K22*K42*r_adab;
	float P__25 = P_25-g0*(K20*P_05+K21*P_15+K22*P_25+K50*P_20+K51*P_21+K52*P_22)+g0_2*(K50*(K20*P_00+K21*P_10+K22*P_20)+K51*(K20*P_01+K21*P_11+K22*P_21)+K52*(K20*P_02+K21*P_12+K22*P_22))+K20*K50*r_adab+K21*K51*r_adab+K22*K52*r_adab;
	float P__30 = P_30-g0*(K00*P_30+K01*P_31+K02*P_32+K30*P_00+K31*P_10+K32*P_20)+g0_2*(K00*(K30*P_00+K31*P_10+K32*P_20)+K01*(K30*P_01+K31*P_11+K32*P_21)+K02*(K30*P_02+K31*P_12+K32*P_22))+K00*K30*r_adab+K01*K31*r_adab+K02*K32*r_adab;
	float P__31 = P_31-g0*(K10*P_30+K11*P_31+K12*P_32+K30*P_01+K31*P_11+K32*P_21)+g0_2*(K10*(K30*P_00+K31*P_10+K32*P_20)+K11*(K30*P_01+K31*P_11+K32*P_21)+K12*(K30*P_02+K31*P_12+K32*P_22))+K10*K30*r_adab+K11*K31*r_adab+K12*K32*r_adab;
	float P__32 = P_32-g0*(K20*P_30+K21*P_31+K22*P_32+K30*P_02+K31*P_12+K32*P_22)+g0_2*(K20*(K30*P_00+K31*P_10+K32*P_20)+K21*(K30*P_01+K31*P_11+K32*P_21)+K22*(K30*P_02+K31*P_12+K32*P_22))+K20*K30*r_adab+K21*K31*r_adab+K22*K32*r_adab;
	float P__33 = P_33-g0*(K30*P_03+K31*P_13+K30*P_30+K31*P_31+K32*P_23+K32*P_32)+(K30*K30)*r_adab+(K31*K31)*r_adab+(K32*K32)*r_adab+g0_2*(K30*(K30*P_00+K31*P_10+K32*P_20)+K31*(K30*P_01+K31*P_11+K32*P_21)+K32*(K30*P_02+K31*P_12+K32*P_22));
	float P__34 = P_34-g0*(K30*P_04+K31*P_14+K32*P_24+K40*P_30+K41*P_31+K42*P_32)+g0_2*(K40*(K30*P_00+K31*P_10+K32*P_20)+K41*(K30*P_01+K31*P_11+K32*P_21)+K42*(K30*P_02+K31*P_12+K32*P_22))+K30*K40*r_adab+K31*K41*r_adab+K32*K42*r_adab;
	float P__35 = P_35-g0*(K30*P_05+K31*P_15+K32*P_25+K50*P_30+K51*P_31+K52*P_32)+g0_2*(K50*(K30*P_00+K31*P_10+K32*P_20)+K51*(K30*P_01+K31*P_11+K32*P_21)+K52*(K30*P_02+K31*P_12+K32*P_22))+K30*K50*r_adab+K31*K51*r_adab+K32*K52*r_adab;
	float P__40 = P_40-g0*(K00*P_40+K01*P_41+K02*P_42+K40*P_00+K41*P_10+K42*P_20)+g0_2*(K00*(K40*P_00+K41*P_10+K42*P_20)+K01*(K40*P_01+K41*P_11+K42*P_21)+K02*(K40*P_02+K41*P_12+K42*P_22))+K00*K40*r_adab+K01*K41*r_adab+K02*K42*r_adab;
	float P__41 = P_41-g0*(K10*P_40+K11*P_41+K12*P_42+K40*P_01+K41*P_11+K42*P_21)+g0_2*(K10*(K40*P_00+K41*P_10+K42*P_20)+K11*(K40*P_01+K41*P_11+K42*P_21)+K12*(K40*P_02+K41*P_12+K42*P_22))+K10*K40*r_adab+K11*K41*r_adab+K12*K42*r_adab;
	float P__42 = P_42-g0*(K20*P_40+K21*P_41+K22*P_42+K40*P_02+K41*P_12+K42*P_22)+g0_2*(K20*(K40*P_00+K41*P_10+K42*P_20)+K21*(K40*P_01+K41*P_11+K42*P_21)+K22*(K40*P_02+K41*P_12+K42*P_22))+K20*K40*r_adab+K21*K41*r_adab+K22*K42*r_adab;
	float P__43 = P_43-g0*(K30*P_40+K31*P_41+K32*P_42+K40*P_03+K41*P_13+K42*P_23)+g0_2*(K30*(K40*P_00+K41*P_10+K42*P_20)+K31*(K40*P_01+K41*P_11+K42*P_21)+K32*(K40*P_02+K41*P_12+K42*P_22))+K30*K40*r_adab+K31*K41*r_adab+K32*K42*r_adab;
	float P__44 = P_44-g0*(K40*P_04+K41*P_14+K40*P_40+K42*P_24+K41*P_41+K42*P_42)+(K40*K40)*r_adab+(K41*K41)*r_adab+(K42*K42)*r_adab+g0_2*(K40*(K40*P_00+K41*P_10+K42*P_20)+K41*(K40*P_01+K41*P_11+K42*P_21)+K42*(K40*P_02+K41*P_12+K42*P_22));
	float P__45 = P_45-g0*(K40*P_05+K41*P_15+K42*P_25+K50*P_40+K51*P_41+K52*P_42)+g0_2*(K50*(K40*P_00+K41*P_10+K42*P_20)+K51*(K40*P_01+K41*P_11+K42*P_21)+K52*(K40*P_02+K41*P_12+K42*P_22))+K40*K50*r_adab+K41*K51*r_adab+K42*K52*r_adab;
	float P__50 = P_50-g0*(K00*P_50+K01*P_51+K02*P_52+K50*P_00+K51*P_10+K52*P_20)+g0_2*(K00*(K50*P_00+K51*P_10+K52*P_20)+K01*(K50*P_01+K51*P_11+K52*P_21)+K02*(K50*P_02+K51*P_12+K52*P_22))+K00*K50*r_adab+K01*K51*r_adab+K02*K52*r_adab;
	float P__51 = P_51-g0*(K10*P_50+K11*P_51+K12*P_52+K50*P_01+K51*P_11+K52*P_21)+g0_2*(K10*(K50*P_00+K51*P_10+K52*P_20)+K11*(K50*P_01+K51*P_11+K52*P_21)+K12*(K50*P_02+K51*P_12+K52*P_22))+K10*K50*r_adab+K11*K51*r_adab+K12*K52*r_adab;
	float P__52 = P_52-g0*(K20*P_50+K21*P_51+K22*P_52+K50*P_02+K51*P_12+K52*P_22)+g0_2*(K20*(K50*P_00+K51*P_10+K52*P_20)+K21*(K50*P_01+K51*P_11+K52*P_21)+K22*(K50*P_02+K51*P_12+K52*P_22))+K20*K50*r_adab+K21*K51*r_adab+K22*K52*r_adab;
	float P__53 = P_53-g0*(K30*P_50+K31*P_51+K32*P_52+K50*P_03+K51*P_13+K52*P_23)+g0_2*(K30*(K50*P_00+K51*P_10+K52*P_20)+K31*(K50*P_01+K51*P_11+K52*P_21)+K32*(K50*P_02+K51*P_12+K52*P_22))+K30*K50*r_adab+K31*K51*r_adab+K32*K52*r_adab;
	float P__54 = P_54-g0*(K40*P_50+K41*P_51+K42*P_52+K50*P_04+K51*P_14+K52*P_24)+g0_2*(K40*(K50*P_00+K51*P_10+K52*P_20)+K41*(K50*P_01+K51*P_11+K52*P_21)+K42*(K50*P_02+K51*P_12+K52*P_22))+K40*K50*r_adab+K41*K51*r_adab+K42*K52*r_adab;
	float P__55 = P_55-g0*(K50*P_05+K51*P_15+K52*P_25+K50*P_50+K51*P_51+K52*P_52)+(K50*K50)*r_adab+(K51*K51)*r_adab+(K52*K52)*r_adab+g0_2*(K50*(K50*P_00+K51*P_10+K52*P_20)+K51*(K50*P_01+K51*P_11+K52*P_21)+K52*(K50*P_02+K51*P_12+K52*P_22));


	float len = sqrt(x0*x0+x1*x1+x2*x2);
	float invlen3 = 1.0/(len*len*len);
	float invlen32 = (invlen3*invlen3);

	float x1_x2 = (x1*x1+x2*x2);
	float x0_x2 = (x0*x0+x2*x2);
	float x0_x1 = (x0*x0+x1*x1);

	// normalized a posteriori covariance
	P00 = invlen32*(-x1_x2*(-P__00*x1_x2+P__10*x0*x1+P__20*x0*x2)+x0*x1*(-P__01*x1_x2+P__11*x0*x1+P__21*x0*x2)+x0*x2*(-P__02*x1_x2+P__12*x0*x1+P__22*x0*x2));
	P01 = invlen32*(-x0_x2*(-P__01*x1_x2+P__11*x0*x1+P__21*x0*x2)+x0*x1*(-P__00*x1_x2+P__10*x0*x1+P__20*x0*x2)+x1*x2*(-P__02*x1_x2+P__12*x0*x1+P__22*x0*x2));
	P02 = invlen32*(-x0_x1*(-P__02*x1_x2+P__12*x0*x1+P__22*x0*x2)+x0*x2*(-P__00*x1_x2+P__10*x0*x1+P__20*x0*x2)+x1*x2*(-P__01*x1_x2+P__11*x0*x1+P__21*x0*x2));
	P03 = -invlen3*(-P__03*x1_x2+P__13*x0*x1+P__23*x0*x2);
	P04 = -invlen3*(-P__04*x1_x2+P__14*x0*x1+P__24*x0*x2);
	P05 = -invlen3*(-P__05*x1_x2+P__15*x0*x1+P__25*x0*x2);
	P10 = invlen32*(-x1_x2*(-P__10*x0_x2+P__00*x0*x1+P__20*x1*x2)+x0*x1*(-P__11*x0_x2+P__01*x0*x1+P__21*x1*x2)+x0*x2*(-P__12*x0_x2+P__02*x0*x1+P__22*x1*x2));
	P11 = invlen32*(-x0_x2*(-P__11*x0_x2+P__01*x0*x1+P__21*x1*x2)+x0*x1*(-P__10*x0_x2+P__00*x0*x1+P__20*x1*x2)+x1*x2*(-P__12*x0_x2+P__02*x0*x1+P__22*x1*x2));
	P12 = invlen32*(-x0_x1*(-P__12*x0_x2+P__02*x0*x1+P__22*x1*x2)+x0*x2*(-P__10*x0_x2+P__00*x0*x1+P__20*x1*x2)+x1*x2*(-P__11*x0_x2+P__01*x0*x1+P__21*x1*x2));
	P13 = -invlen3*(-P__13*x0_x2+P__03*x0*x1+P__23*x1*x2);
	P14 = -invlen3*(-P__14*x0_x2+P__04*x0*x1+P__24*x1*x2);
	P15 = -invlen3*(-P__15*x0_x2+P__05*x0*x1+P__25*x1*x2);
	P20 = invlen32*(-x1_x2*(-P__20*x0_x1+P__00*x0*x2+P__10*x1*x2)+x0*x1*(-P__21*x0_x1+P__01*x0*x2+P__11*x1*x2)+x0*x2*(-P__22*x0_x1+P__02*x0*x2+P__12*x1*x2));
	P21 = invlen32*(-x0_x2*(-P__21*x0_x1+P__01*x0*x2+P__11*x1*x2)+x0*x1*(-P__20*x0_x1+P__00*x0*x2+P__10*x1*x2)+x1*x2*(-P__22*x0_x1+P__02*x0*x2+P__12*x1*x2));
	P22 = invlen32*(-x0_x1*(-P__22*x0_x1+P__02*x0*x2+P__12*x1*x2)+x0*x2*(-P__20*x0_x1+P__00*x0*x2+P__10*x1*x2)+x1*x2*(-P__21*x0_x1+P__01*x0*x2+P__11*x1*x2));
	P23 = -invlen3*(-P__23*x0_x1+P__03*x0*x2+P__13*x1*x2);
	P24 = -invlen3*(-P__24*x0_x1+P__04*x0*x2+P__14*x1*x2);
	P25 = -invlen3*(-P__25*x0_x1+P__05*x0*x2+P__15*x1*x2);
	P30 = -invlen3*(-P__30*x1_x2+P__31*x0*x1+P__32*x0*x2);
	P31 = -invlen3*(-P__31*x0_x2+P__30*x0*x1+P__32*x1*x2);
	P32 = -invlen3*(-P__32*x0_x1+P__30*x0*x2+P__31*x1*x2);
	P33 = P__33;
	P34 = P__34;
	P35 = P__35;
	P40 = -invlen3*(-P__40*x1_x2+P__41*x0*x1+P__42*x0*x2);
	P41 = -invlen3*(-P__41*x0_x2+P__40*x0*x1+P__42*x1*x2);
	P42 = -invlen3*(-P__42*x0_x1+P__40*x0*x2+P__41*x1*x2);
	P43 = P__43;
	P44 = P__44;
	P45 = P__45;
	P50 = -invlen3*(-P__50*x1_x2+P__51*x0*x1+P__52*x0*x2);
	P51 = -invlen3*(-P__51*x0_x2+P__50*x0*x1+P__52*x1*x2);
	P52 = -invlen3*(-P__52*x0_x1+P__50*x0*x2+P__51*x1*x2);
	P53 = P__53;
	P54 = P__54;
	P55 = P__55;

	// normalized a posteriori state
	x0 = x0/len;
	x1 = x1/len;
	x2 = x2/len;


	// compute Euler angles (not exactly a part of the extended Kalman filter)
	// yaw integration through full rotation matrix
	float u_nb1 = u1 - x4;
	float u_nb2 = u2 - x5;

	float cy = cos(yaw); //old angles (last state before integration)
	float sy = sin(yaw);
	float d = sqrt(x_last[1]*x_last[1] + x_last[2]*x_last[2]);
	float d_inv = 1.0f / d;

	// compute needed parts of rotation matrix R (state and angle based version, equivalent with the commented version above)
	float R11 = cy * d;
	float R12 = -(x_last[2]*sy + x_last[0]*x_last[1]*cy) * d_inv;
	float R13 = (x_last[1]*sy - x_last[0]*x_last[2]*cy) * d_inv;
	float R21 = sy * d;
	float R22 = (x_last[2]*cy - x_last[0]*x_last[1]*sy) * d_inv;
	float R23 = -(x_last[1]*cy + x_last[0]*x_last[2]*sy) * d_inv;

	// update needed parts of R for yaw computation
	float R11_new = R11 + h*(u_nb2*R12 - u_nb1*R13);
	float R21_new = R21 + h*(u_nb2*R22 - u_nb1*R23);
	yaw = atan2(R21_new,R11_new);

	// compute new pitch and roll angles from a posteriori states
	pitch = asin(-x0);
	roll = atan2(x1,x2);

	// save the estimated non-gravitational acceleration
	a0 = z0-x0*g0;
	a1 = z1-x1*g0;
	a2 = z2-x2*g0;
}

void DCM_IMU_uC::getState(float *State) {
	State[0] = x0;
	State[1] = x1;
	State[2] = x2;
	State[3] = x3;
	State[4] = x4;
	State[5] = x5;
}

void DCM_IMU_uC::getCovariance(float *Covariance) {
	Covariance[0] = P00;
	Covariance[1] = P01;
	Covariance[2] = P02;
	Covariance[3] = P03;
	Covariance[4] = P04;
	Covariance[5] = P05;
	Covariance[6] = P10;
	Covariance[7] = P11;
	Covariance[8] = P12;
	Covariance[9] = P13;
	Covariance[10] = P14;
	Covariance[11] = P15;
	Covariance[12] = P20;
	Covariance[13] = P21;
	Covariance[14] = P22;
	Covariance[15] = P23;
	Covariance[16] = P24;
	Covariance[17] = P25;
	Covariance[18] = P30;
	Covariance[19] = P31;
	Covariance[20] = P32;
	Covariance[21] = P33;
	Covariance[22] = P34;
	Covariance[23] = P35;
	Covariance[24] = P40;
	Covariance[25] = P41;
	Covariance[26] = P42;
	Covariance[27] = P43;
	Covariance[28] = P44;
	Covariance[29] = P45;
	Covariance[30] = P50;
	Covariance[31] = P51;
	Covariance[32] = P52;
	Covariance[33] = P53;
	Covariance[34] = P54;
	Covariance[35] = P55;
}

void DCM_IMU_uC::getNGAcc(float *ngacc) {
	ngacc[0] = a0;
	ngacc[1] = a1;
	ngacc[2] = a2;
}

