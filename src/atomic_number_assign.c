#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int AtomicNumber(char data[3])
{
	struct Tabla
	{	
		char *name_element;
        	int atomic_number;
	};

	struct Tabla element[97]; 

	//******** Element's name**********************
	element[0].name_element = "";
	element[1].name_element = "H";
	element[2].name_element = "He";
	element[3].name_element = "Li";
	element[4].name_element = "Be";
	element[5].name_element = "B";
	element[6].name_element = "C";
	element[7].name_element = "N";
	element[8].name_element = "O";
	element[9].name_element = "F";
	element[10].name_element = "Ne";
	element[11].name_element = "Na";
	element[12].name_element = "Mg";
	element[13].name_element = "Al";
	element[14].name_element = "Si";
	element[15].name_element = "P";
	element[16].name_element = "S";
	element[17].name_element = "Cl";
	element[18].name_element = "Ar";
	element[19].name_element = "K";
	element[20].name_element = "Ca";
	element[21].name_element = "Sc";
	element[22].name_element = "Ti";
	element[23].name_element = "V";
	element[24].name_element = "Cr";
	element[25].name_element = "Mn";
	element[26].name_element = "Fe";
	element[27].name_element = "Co";
	element[28].name_element = "Ni";
	element[29].name_element = "Cu";
	element[30].name_element = "Zn";
	element[31].name_element = "Ga";
	element[32].name_element = "Ge";
	element[33].name_element = "As";
	element[34].name_element = "Se";
	element[35].name_element = "Br";
	element[36].name_element = "Kr";
	element[37].name_element = "Rb";
	element[38].name_element = "Sr";
	element[39].name_element = "Y";
	element[40].name_element = "Zr";
	element[41].name_element = "Nb";
	element[42].name_element = "Mo";
	element[43].name_element = "Tc";
	element[44].name_element = "Ru";
	element[45].name_element = "Rh";
	element[46].name_element = "Pd";
	element[47].name_element = "Ag";
	element[48].name_element = "Cd";
	element[49].name_element = "In";
	element[50].name_element = "Sn";
	element[51].name_element = "Sb";
	element[52].name_element = "Te";
	element[53].name_element = "I";
	element[54].name_element = "Xe";
	element[55].name_element = "Cs";
	element[56].name_element = "Ba";
	element[57].name_element = "La";
	element[58].name_element = "Ce";
	element[59].name_element = "Pr";
	element[60].name_element = "Nd";
	element[61].name_element = "Pm";
	element[62].name_element = "Sm";
	element[63].name_element = "Eu";
	element[64].name_element = "Gd";
	element[65].name_element = "Tb";
	element[66].name_element = "Dy";
	element[67].name_element = "Ho";
	element[68].name_element = "Er";
	element[69].name_element = "Tm";
	element[70].name_element = "Yb";
	element[71].name_element = "Lu";
	element[72].name_element = "Hf";
	element[73].name_element = "Ta";
	element[74].name_element = "W";
	element[75].name_element = "Re";
	element[76].name_element = "Os";
	element[77].name_element = "Ir";
	element[78].name_element = "Pt";
	element[79].name_element = "Au";
	element[80].name_element = "Hg";
	element[81].name_element = "Tl";
	element[82].name_element = "Pb";
	element[83].name_element = "Bi";
	element[84].name_element = "Po";
	element[85].name_element = "At";
	element[86].name_element = "Rn";
	element[87].name_element = "Fr";
	element[88].name_element = "Ra";
	element[89].name_element = "Ac";
	element[90].name_element = "Th";
	element[91].name_element = "Pa";
	element[92].name_element = "U";
	element[93].name_element = "Np";
	element[94].name_element = "Pu";
	element[95].name_element = "Am";
	element[96].name_element = "Cm";
	//******** End of Element's name************************
	//
	//
	//
	//******** Element's Atomic Number **********************
	
	element[0].atomic_number = 0;
	element[1].atomic_number = 1;
	element[2].atomic_number = 2;
	element[3].atomic_number = 3;
	element[4].atomic_number = 4;
	element[5].atomic_number = 5;
	element[6].atomic_number = 6;
	element[7].atomic_number = 7;
	element[8].atomic_number = 8;
	element[9].atomic_number = 9;
	element[10].atomic_number = 10;
	element[11].atomic_number = 11;
	element[12].atomic_number = 12;
	element[13].atomic_number = 13;
	element[14].atomic_number = 14;
	element[15].atomic_number = 15;
	element[16].atomic_number = 16;
	element[17].atomic_number = 17;
	element[18].atomic_number = 18;
	element[19].atomic_number = 19;
	element[20].atomic_number = 20;
	element[21].atomic_number = 21;
	element[22].atomic_number = 22;
	element[23].atomic_number = 23;
	element[24].atomic_number = 24;
	element[25].atomic_number = 25;
	element[26].atomic_number = 26;
	element[27].atomic_number = 27;
	element[28].atomic_number = 28;
	element[29].atomic_number = 29;
	element[30].atomic_number = 30;
	element[31].atomic_number = 31;
	element[32].atomic_number = 32;
	element[33].atomic_number = 33;
	element[34].atomic_number = 34;
	element[35].atomic_number = 35;
	element[36].atomic_number = 36;
	element[37].atomic_number = 37;
	element[38].atomic_number = 38;
	element[39].atomic_number = 39;
	element[40].atomic_number = 40;
	element[41].atomic_number = 41;
	element[42].atomic_number = 42;
	element[43].atomic_number = 43;
	element[44].atomic_number = 44;
	element[45].atomic_number = 45;
	element[46].atomic_number = 46;
	element[47].atomic_number = 47;
	element[48].atomic_number = 48;
	element[49].atomic_number = 49;
	element[50].atomic_number = 50;
	element[51].atomic_number = 51;
	element[52].atomic_number = 52;
	element[53].atomic_number = 53;
	element[54].atomic_number = 54;
	element[55].atomic_number = 55;
	element[56].atomic_number = 56;
	element[57].atomic_number = 57;
	element[58].atomic_number = 58;
	element[59].atomic_number = 59;
	element[60].atomic_number = 60;
	element[61].atomic_number = 61;
	element[62].atomic_number = 62;
	element[63].atomic_number = 63;
	element[64].atomic_number = 64;
	element[65].atomic_number = 65;
	element[66].atomic_number = 66;
	element[67].atomic_number = 67;
	element[68].atomic_number = 68;
	element[69].atomic_number = 69;
	element[70].atomic_number = 70;
	element[71].atomic_number = 71;
	element[72].atomic_number = 72;
	element[73].atomic_number = 73;
	element[74].atomic_number = 74;
	element[75].atomic_number = 75;
	element[76].atomic_number = 76;
	element[77].atomic_number = 77;
	element[78].atomic_number = 78;
	element[79].atomic_number = 79;
	element[80].atomic_number = 80;
	element[81].atomic_number = 81;
	element[82].atomic_number = 82;
	element[83].atomic_number = 83;
	element[84].atomic_number = 84;
	element[85].atomic_number = 85;
	element[86].atomic_number = 86;
	element[87].atomic_number = 87;
	element[88].atomic_number = 88;
	element[89].atomic_number = 89;
	element[90].atomic_number = 90;
	element[91].atomic_number = 91;
	element[92].atomic_number = 92;
	element[93].atomic_number = 93;
	element[94].atomic_number = 94;
	element[95].atomic_number = 95;
	element[96].atomic_number = 96;
	//******** End of Element's Atomic Number************************

	int target;
	for (int i=0; i<97; i++)
	{
		if(strcmp(data,element[i].name_element)==0)
		{
			target = element[i].atomic_number;
			break;
		}
	}
	return target;

}
