#include "wla.h"
#include <iostream>//Just for testing purposes...
#include <cmath>

double awc(int, double *);
void linear_regression(int, double *, double *);

dwt::dwt(int width, int type){
  w = width;
  family = type;
  g = new double[w];
  h = new double[w];
  init_coeff(w, family, g, h);
}

dwt::~dwt(){
  delete[] g;
  delete[] h;
}

void dwt::init_coeff(int w, int family, double *g, double *h){
  if(family == 1){              //Make switch case instead of if else?
    for(int i = 0; i < w; i++){
      g[i] = get_daub(w, i);
    }
  }else if(family == 2){
    for(int i = 0; i < w; i++){
      g[i] = get_la(w, i);
    }
  }else if(family == 3){
    for(int i = 0; i < w; i++){
      g[i] = get_bl(w, i);
    }
  }else if(family == 4){
    for(int i = 0; i < w; i++){
      g[i] = get_coiflet(w, i);
    }
  }
  
  for(int i = 0; i < w; i++){
    h[i] = pow(-1,i)*g[w-1-i];
  }
}

void dwt::transform(int n, double *v_in, double *v_out, double *w_out){
  int k;
  int nh = n/2;

  for(int i = 0; i < nh; i++){
    k = 2*i+1;
    v_out[i] = g[0]*v_in[k];
    w_out[i] = h[0]*v_in[k];
    for(int j = 1; j < w; j++){
      k--;
      if(k < 0){
	k = n-1;
      }
      v_out[i] += g[j]*v_in[k];
      w_out[i] += h[j]*v_in[k];
    }
  }
}

void dwt::inverse_transform(int n, double *v_in, double *w_in, double *x_out){
  int k, l, m;
  int o = -2; 
  int p = -1;  
  int wh = w/2;
 
  for(int i = 0; i < n; i++){
    k = i;
    l = 0;
    m = 1;
    o += 2;
    p += 2;
    x_out[o] = h[m]*w_in[k]+g[m]*v_in[k];
    x_out[p] = h[l]*w_in[k]+g[l]*v_in[k];
    if(w > 2){
      for(int j = 1; j < wh; j++){
	k++;
	if(k >= n){
	  k = 0;
	}
	l += 2;
	m += 2;
	x_out[o] += h[m]*w_in[k]+g[m]*v_in[k];
	x_out[p] += h[l]*w_in[k]+g[l]*v_in[k];
      }
    }
  }
}

modwt::modwt(int width, int type){
  w = width;
  family = type;
  it = 0;//Create another constructor which sets this value to 0?
  g = new double[w];
  h = new double[w];
  init_coeff(w, family, g, h);
}

modwt::~modwt(){
  delete[] g;
  delete[] h;
}

void modwt::init_coeff(int w, int family, double *g, double *h){
  if(family == 1){              //Make switch case instead of if else?
    for(int i = 0; i < w; i++){
      g[i] = get_daub(w, i)/sqrt(2);
    }
  }else if(family == 2){
    for(int i = 0; i < w; i++){
      g[i] = get_la(w, i)/sqrt(2);
    }
  }else if(family == 3){
    for(int i = 0; i < w; i++){
      g[i] = get_bl(w, i)/sqrt(2);
    }
  }else if(family == 4){
    for(int i = 0; i < w; i++){
      g[i] = get_coiflet(w, i)/sqrt(2);
    }
  }
  
  for(int i = 0; i < w; i++){
    h[i] = pow(-1,i)*g[w-1-i];
  }
}

void modwt::transform(int n, double *v_in, double *v_out, double *w_out){
  int k;
  
  it++;                       //iteration... Watch out if I run this code multiple times  
  for(int i = 0; i < n; i++){
    k = i;
    v_out[i] = g[0]*v_in[k];
    w_out[i] = h[0]*v_in[k];
    for(int j = 1; j < w; j++){
      k -= pow(2,it-1);
      if(k < 0){
	k += n;
      }
      v_out[i] += g[j]*v_in[k];
      w_out[i] += h[j]*v_in[k];
    }
  }
}

void modwt::inverse_transform(int n, double *v_in, double *w_in, double *x_out){
  int k;
  
  for(int i = 0; i < n; i++){
    k = i;
    x_out[i] = h[0]*w_in[k]+g[0]*v_in[k];
    for(int j = 1; j < w; j++){
      k += pow(2, it-1);
      if(k >= n){
	k -= n;
      }
      x_out[i] += h[j]*w_in[k]+g[j]*v_in[k];
    }
  }
  it--; 
}

lifting::lifting(){
  //Constructor for the lifting algorithm
}

cwt::cwt(){
  //Constructor for the cwt
}

double get_daub(int l, int i){ //Er dette dumt?
  const double d2[2] = {            //Lager ikke programmet nye arrays her hver gang denne funksjonen kalles
    0.7071067811865475,             //Gj�re dem globale istedet??? static???
    0.7071067811865475
  };
  const double d4[4] = {
    0.4829629131445341,
    0.8365163037378077,
    0.2241438680420134,
    -0.1294095225512603
  };
  const double d6[6] = {
    0.3326705529500827,
    0.8068915093110928,
    0.4598775021184915,
    -0.1350110200102546,
    -0.0854412738820267,
    0.0352262918857096,
  };
  const double d8[8] = {
    0.2303778133074431,
    0.7148465705484058,
    0.6308807679358788,
    -0.0279837694166834,
    -0.1870348117179132,
    0.0308413818353661,
    0.0328830116666778,
    -0.0105974017850021
  };
  const double d10[10] = {
    0.1601023979741930,
    0.6038292697971898,
    0.7243085284377729,
    0.1384281459013204,
    -0.2422948870663824,
    -0.0322448695846381,
    0.0775714938400459,
    -0.0062414902127983,
    -0.0125807519990820,
    0.0033357252854738
  };
  const double d12[12] = {
    0.1115407433501094,
    0.4946238903984530,
    0.7511339080210954,
    0.3152503517091980,
    -0.2262646939654399,
    -0.1297668675672624,
    0.0975016055873224,
    0.0275228655303053,
    -0.0315820393174862,
    0.0005538422011614,
    0.0047772575109455,
    -0.0010773010853085
  };
  const double d14[14] = {
    0.0778520540850081,
    0.3965393194819136,
    0.7291320908462368,
    0.4697822874052154,
    -0.1439060039285293,
    -0.2240361849938538,
    0.0713092192668312,
    0.0806126091510820,
    -0.0380299369350125,
    -0.0165745416306664,
    0.0125509985560993,
    0.0004295779729214,
    -0.0018016407040474,
    0.0003537137999745
  };
  const double d16[16] = {
    0.0544158422431049,
    0.3128715909143031,
    0.6756307362972904,
    0.5853546836541907,
    -0.0158291052563816,
    -0.2840155429615702,
    0.0004724845739124,
    0.1287474266204837,
    -0.0173693010018083,
    -0.0440882539307952,
    0.0139810279173995,
    0.0087460940474061,
    -0.0048703529934518,
    -0.0003917403733770,
    0.0006754494064506,
    -0.0001174767841248
  };
  const double d18[18] = {
    0.0380779473638791,
    0.2438346746125939,
    0.6048231236901156,
    0.6572880780512955,
    0.1331973858249927,
    -0.2932737832791761,
    -0.0968407832229524,
    0.1485407493381306,
    0.0307256814793395,
    -0.0676328290613302,
    0.0002509471148340,
    0.0223616621236805,
    -0.0047232047577520,
    -0.0042815036824636,
    0.0018476468830564,
    0.0002303857635232,
    -0.0002519631889427,
    0.0000393473203163
  };
  const double d20[20] = {
    0.0266700579005546,
    0.1881768000776863,
    0.5272011889317202,
    0.6884590394536250,
    0.2811723436606485,
    -0.2498464243272283,
    -0.1959462743773399,
    0.1273693403357890,
    0.0930573646035802,
    -0.0713941471663697,
    -0.0294575368218480,
    0.0332126740593703,
    0.0036065535669880,
    -0.0107331754833036,
    0.0013953517470692,
    0.0019924052951930,
    -0.0006858566949566,
    -0.0001164668551285,
    0.0000935886703202,
    -0.0000132642028945
    };

  if(l == 2){
    return d2[i];
  }else if(l == 4){
    return d4[i];
  }else if(l == 6){
    return d6[i];
  }else if(l == 8){
    return d8[i];
  }else if(l == 10){
    return d10[i];
  }else if(l == 12){
    return d12[i];
  }else if(l == 14){
    return d14[i];
  }else if(l == 16){
    return d16[i];
  }else if(l == 18){
    return d18[i];
  }else if(l == 20){
    return d20[i];
  }
}

double get_la(int l, int i){
  const double la8[8] = {
    -0.0757657147893407,
    -0.0296355276459541,
    0.4976186676324578,
    0.8037387518052163,
    0.2978577956055422,
    -0.0992195435769354,
    -0.0126039672622612,
    0.0322231006040713
  };
  const double la10[10] = {
    0.0195388827353869,
    -0.0211018340249298,
    -0.1753280899081075,
    0.0166021057644243,
    0.6339789634569490,
    0.7234076904038076,
    0.1993975339769955,
    -0.0391342493025834,
    0.0295194909260734,
    0.0273330683451645
  };
  const double la12[12] = {
    0.0154041093273377,
    0.0034907120843304,
    -0.1179901111484105,
    -0.0483117425859981,
    0.4910559419276396,
    0.7876411410287941,
    0.3379294217282401,
    -0.0726375227866000,
    -0.0210602925126954,
    0.0447249017707482,
    0.0017677118643983,
    -0.0078007083247650
  };
  const double la14[14] = {
    0.0102681767084968,
    0.0040102448717033,
    -0.1078082377036168,
    -0.1400472404427030,
    0.2886296317509833,
    0.7677643170045710,
    0.5361019170907720,
    0.0174412550871099,
    -0.0495528349370410,
    0.0678926935015971,
    0.0305155131659062,
    -0.0126363034031526,
    -0.0010473848889657,
    0.0026818145681164
  };
  const double la16[16] = {
    -0.0033824159513594,
    -0.0005421323316355,
    0.0316950878103452,
    0.0076074873252848,
    -0.1432942383510542,
    -0.0612733590679088,
    0.4813596512592012,
    0.7771857516997478,
    0.3644418948359564,
    -0.0519458381078751,
    -0.0272190299168137,
    0.0491371796734768,
    0.0038087520140601,
    -0.0149522583367926,
    -0.0003029205145516,
    0.0018899503329007
  };
  const double la18[18] = {
    0.0010694900326538,
    -0.0004731544985879,
    -0.0102640640276849,
    0.0088592674935117,
    0.0620777893027638,
    -0.0182337707798257,
    -0.1915508312964873,
    0.0352724880359345,
    0.6173384491413523,
    0.7178970827642257,
    0.2387609146074182,
    -0.0545689584305765,
    0.0005834627463312,
    0.0302248788579895,
    -0.0115282102079848,
    -0.0132719677815332,
    0.0006197808890549,
    0.0014009155255716
  };
  const double la20[20] = {
    0.0007701598091030,
    0.0000956326707837,
    -0.0086412992759401,
    -0.0014653825833465,
    0.0459272392237649,
    0.0116098939129724,
    -0.1594942788575307,
    -0.0708805358108615,
    0.4716906668426588,
    0.7695100370143388,
    0.3838267612253823,
    -0.0355367403054689,
    -0.0319900568281631,
    0.0499949720791560,
    0.0057649120455518,
    -0.0203549398039460,
    -0.0008043589345370,
    0.0045931735836703,
    0.0000570360843390,
    -0.0004593294205481
  };

  if(l == 8){
    return la8[i];
  }else if(l == 10){
    return la10[i];
  }else if(l == 12){
    return la12[i];
  }else if(l == 14){
    return la14[i];
  }else if(l == 16){
    return la16[i];
  }else if(l == 18){
    return la18[i];
  }else if(l == 20){
    return la20[i];
  }
}

double get_bl(int l, int i){
  const double bl14[14] = {
    0.0120154192834842,
    0.0172133762994439,
    -0.0649080035533744,
    -0.0641312898189170,
    0.3602184608985549,
    0.7819215932965554,
    0.4836109156937821,
    -0.0568044768822707,
    -0.1010109208664125,
    0.0447423494687405,
    0.0204642075778225,
    -0.0181266051311065,
    -0.0032832978473081,
    0.0022918339541009
  };
  const double bl18[18] = {
    0.0002594576266544,
    -0.0006273974067728,
    -0.0019161070047557,
    0.0059845525181721,
    0.0040676562965785,
    -0.0295361433733604,
    -0.0002189514157348,
    0.0856124017265279,
    -0.0211480310688774,
    -0.1432929759396520,
    0.2337782900224977,
    0.7374707619933686,
    0.5926551374433956,
    0.0805670008868546,
    -0.1143343069619310,
    -0.0348460237698368,
    0.0139636362487191,
    0.0057746045512475
  };
  const double bl20[20] = {
    0.0008625782242896,
    0.0007154205305517,
    -0.0070567640909701,
    0.0005956827305406,
    0.0496861265075979,
    0.0262403647054251,
    -0.1215521061578162,
    -0.0150192395413644,
    0.5137098728334054,
    0.7669548365010849,
    0.3402160135110789,
    -0.0878787107378667,
    -0.0670899071680668,
    0.0338423550064691,
    -0.0008687519578684,
    -0.0230054612862905,
    -0.0011404297773324,
    0.0050716491945793,
    0.0003401492622332,
    -0.0004101159165852
  };

  if(l == 14){
    return bl14[i];
  }else if(l == 18){
    return bl18[i];
  }else if(l == 20){
    return bl20[i];
  }
}

double get_coiflet(int l, int i){
  const double c6[6] = {
    -0.0156557285289848,
    -0.0727326213410511,
    0.3848648565381134,
    0.8525720416423900,
    0.3378976709511590,
    -0.0727322757411889
  };
  const double c12[12] = {
    -0.0007205494453679,
    -0.0018232088707116,
    0.0056114348194211,
    0.0236801719464464,
    -0.0594344186467388,
    -0.0764885990786692,
    0.4170051844236707,
    0.8127236354493977,
    0.3861100668229939,
    -0.0673725547222826,
    -0.0414649367819558,
    0.0163873364635998
  };
  const double c18[18] = {
    -0.0000345997728362,
    -0.0000709833031381,
    0.0004662169601129,
    0.0011175187708906,
    -0.0025745176887502,
    -0.0090079761366615,
    0.0158805448636158,
    0.0345550275730615,
    -0.0823019271068856,
    -0.0717998216193117,
    0.4284834763776168,
    0.7937772226256169,
    0.4051769024096150,
    -0.0611233900026726,
    -0.0657719112818552,
    0.0234526961418362,
    0.0077825964273254,
    -0.0037935128644910
  };
  const double c24[24] = {
    -0.0000017849850031,
    -0.0000032596802369,
    0.0000312298758654,
    0.0000623390344610,
    -0.0002599745524878,
    -0.0005890207562444,
    0.0012665619292991,
    0.0037514361572790,
    -0.0056582866866115,
    -0.0152117315279485,
    0.0250822618448678,
    0.0393344271233433,
    -0.0962204420340021,
    -0.0666274742634348,
    0.4343860564915321,
    0.7822389309206135,
    0.4153084070304910,
    -0.0560773133167630,
    -0.0812666996808907,
    0.0266823001560570,
    0.0160689439647787,
    -0.0073461663276432,
    -0.0016294920126020,
    0.0008923136685824
  };
  const double c30[30] = {
    -0.0000000951765727,
    -0.0000001674428858,
    0.0000020637618516,
    0.0000037346551755,
    -0.0000213150268122,
    -0.0000413404322768,
    0.0001405411497166,
    0.0003022595818445,
    -0.0006381313431115,
    -0.0016628637021860,
    0.0024333732129107,
    0.0067641854487565,
    -0.0091642311634348,
    -0.0197617789446276,
    0.0326835742705106,
    0.0412892087544753,
    -0.1055742087143175,
    -0.0620359639693546,
    0.4379916262173834,
    0.7742896037334738,
    0.4215662067346898,
    -0.0520431631816557,
    -0.0919200105692549,
    0.0281680289738655,
    0.0234081567882734,
    -0.0101311175209033,
    -0.0041593587818186,
    0.0021782363583355,
    0.0003585896879330,
    -0.0002120808398259
  };

  if(l == 6){
    return c6[i];
  }else if(l == 12){
    return c12[i];
  }else if(l == 18){
    return c18[i];
  }else if(l == 24){
    return c24[i];
  }else if(l == 30){
    return c30[i];
  }
}

double awc(int n, double *w){//Hva med � summe w^2/a^2 istedet????????? 1.2 i time series boka...
  double sum = .0;
  
  for(int i = 0; i < n; i++){//Plukker kun ut waveletkoeffisientene. M� endres om jeg bruker u og v
    sum += fabs(w[i]);
  }
  
  return sum/n;
}

void linear_regression(int n, double *y, double *stat_data){
  //stat_data = {beta_0, beta_1, sigma, r, r_sq, H}         
  double mean_x = .0;
  double mean_y = .0;
  double s_x = .0;
  double s_y = .0;
  double r = .0;
  double sigma = .0;
  double beta_0, beta_1;
  double *x = new double[n];//Mer ryddig om denne kommer inn som et argument???
  double *y_hat = new double[n];

  for(int i = 0; i < n; i++){
    x[n-1-i] = log2(1./pow(2,i+1));//Blir dette logaritmetullet riktig????????
    y[i] = log2(y[i]);//log2 kanskje!!!
  }
  
  for(int i = 0; i < n; i++){
    mean_x += x[i];
    mean_y += y[i];
  }
  mean_x /= n;
  mean_y /= n;
  for(int i = 0; i < n; i++){
    s_x += pow(x[i]-mean_x,2);
    s_y += pow(y[i]-mean_x,2);
  }
  s_x /= n-1;
  s_x = sqrt(s_x);
  s_y /= n-1;
  s_y = sqrt(s_y);
  for(int i = 0; i < n; i++){
    r += (x[i]-mean_x)*(y[i]-mean_y)/(s_x*s_y);
  }
  r /= n-1;
  
  beta_1 = r*s_y/s_x;           //M� endre her for � tilpasse logaritmisk skala?
  beta_0 = mean_y-beta_1*mean_x;

  for(int i = 0; i < n; i++){
    y_hat[i] = beta_0+beta_1*x[i];
    sigma += pow(y[i] - y_hat[i],2);
  }
  sigma /= n-2;
  
  stat_data[0] = beta_0;
  stat_data[1] = beta_1;
  stat_data[2] = sigma;
  stat_data[3] = r;
  stat_data[4] = r*r;
  stat_data[5] = beta_1-.5;

  delete[] x;
  delete[] y_hat;
}



/*
void dwt::transform(int n, double *x){//Transform will execute faster if I use two half-sized arrays (v,w)
  int j = 0;
  int nh = n/2;
  double *x_0 = new double[n+w-2];
  
  for(int i = 0; i < w-2; i++){
    x_0[i+n] = x[i]; //What about the boundary conditions??????????
  }
  for(int i = 0; i < n; i++){
    x_0[i] = x[i];
    x[i] = 0;
  }
  
  for(int i = 0; i < nh; i++){
    for(int k = 0; k < w; k++){
      x[i] += g[k]*x_0[j+k];
      x[i+nh] += h[k]*x_0[j+k];
    }
    j+=2;
  }

  delete[] x_0;
}

void dwt::inverse_transform(int n, double *b){
  int j = 0;
  int l = 1;
  int nh = n/2;
  int wh = w/2-1;
  int kh;
  double *x = new double[n];
  
  for(int i = 0; i < n; i++){
    x[i] = b[i];
    b[i] = 0;
  }
  
  for(int i = 0; i < wh; i++){
    for(int k = 0; k < w; k+=2){
      kh = k/2;
      if(k > j){
	b[j] += g[k]*x[nh-l]+h[k]*x[n-l];
	b[j+1] += g[k+1]*x[nh-l]+h[k+1]*x[n-l];
	l++;
      }else{
	b[j] += g[k]*x[i-kh]+h[k]*x[i+nh-kh];
	b[j+1] += g[k+1]*x[i-kh]+h[k+1]*x[i+nh-kh];
	}
    }
    l = 1;
    j += 2;
  }

  for(int i = wh; i < nh; i++){
    for(int k = 0; k < w; k+=2){
      kh = k/2;
      b[j] += g[k]*x[i-kh]+h[k]*x[i+nh-kh];
      b[j+1] += g[k+1]*x[i-kh]+h[k+1]*x[i+nh-kh];
    }
    j += 2;
  }
  
  delete[] x;
}
*/

