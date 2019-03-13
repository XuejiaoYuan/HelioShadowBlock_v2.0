
#include "./sdbkCalculator/sdbkCalculator.h"
#include "./DataStructure/SunRay.h"
#include "./LSPIA/LSF.h"
#include "./DataStructure/Timer.h"


int main(int argc, char** argv) {

	// ������Ӱ���ڵ�����
	// scene_file ���������ļ�
	// sunray_file ̫�������ļ�
	// options: -a_c ����ÿ��ʱ�������о��ӵ���Ӱ�ڵ����(clipper);
	//			-a_r ����ÿ��ʱ�������о��ӵ���Ӱ�ڵ������raytracing��;
	//			-s_l ����ÿ��ʱ�������о��ӵ���Ӱ�ڵ������LSPIA��
	// outfile_options: -f �洢������Ӱ�ڵ��Ľ��(clipper)
	//					-o ���洢���
	// save_path ����ļ��Ĵ洢��ַ
	if (argc != 6) {
		cout << argc << endl;
		cout << "Usage: [scene_file] [sunray_file] [options] [outfile_options] [save_path]" << endl;
		return -1;
	}


	string scene_filepath = string(argv[1]);
	string sunray_filepath = string(argv[2]);
	string options = string(argv[3]);
	string outfile_options = string(argv[4]);
	string save_path = string(argv[5]);

	SunRay sunray;
	Vector3d sunray_dir = sunray.calcSunRay(sunray_filepath);
	MatrixXd a;

	SolarScene* solar_scene = new SolarScene();
	solar_scene->initFieldParam(scene_filepath);

	// Gauss-Legendre����Ԥ����
	int N = solar_scene->recvs[0]->recv_size.x() / RECEIVER_SLICE;
	int M = solar_scene->recvs[0]->recv_size.z() / RECEIVER_SLICE;
	GaussLegendre* gl = new GaussLegendre();
	N = 8;
	M = 8;
	gl->calcNodeWeight(N, M);

	// TODO ���������Ż�
	vector<vector<double>*> field_args;
	switch (solar_scene->layouts[0]->layout_type)
	{
	case CrossRectFieldType:
	case RectLayoutType:
		field_args.push_back(new vector<double>{10, 10, 10});			// ���վ����
		field_args.push_back(new vector<double>{ -80 });					// ��һ�ж��վ��������֮��ľ���
		field_args.push_back(new vector<double>{ 300 });					// ���վ�����
		field_args.push_back(new vector<double>{ 300 });					// ���վ�����
		break;
	case FermatLayoutType:
		field_args.push_back(new vector<double>{ double(0.1) });			// ���վ���Χ�а�ȫ����
		field_args.push_back(new vector<double>{ double( 114.9417 )});	// ��һ��ͬ��Բ�������֮��ľ���(Campo: 114.9417)
		field_args.push_back(new vector<double>{ 0.866f });				// ��һ��ͬ�Ļ��ж��վ��ķֲ����
		field_args.push_back(new vector<double>{ 1.4f });				// �ڶ���ͬ�Ļ��ж��վ��ķֲ����
		field_args.push_back(new vector<double>{ 2.0f});					// ������ͬ�Ļ��ж��վ��ķֲ����
		break;
	default:
		break;
	}
	solar_scene->adjustFieldParam(field_args);

	SdBkCalcCreator sdbk_calc_creator;
	SdBkCalc* sdbk_calc = sdbk_calc_creator.getSdBkCalc(solar_scene, gl);


	// ��С�����������
	//vector<vector<int>> daily = {
	//	{ 3,21 },		// ����
	//	{ 6,22 },		// ����
	//	{ 9,23 },		// ���
	//	{ 12,22 }		// ����
	//};
	//vector<int> sample_hour = { 8, 10, 12, 14, 16 };


	fstream outFile(save_path + "/time_res.txt", ios_base::out);
	vector<int> time_param(4);
	for (int month = 1; month <= 12; month++) {
		for (int day = 1; day < 29; day += 3) {
			for (int hour = 8; hour < 17; hour++) {
				for (int min = 0; min < 30; min += 10) {
					cout << "Calc start" << endl;
					time_param[0] = month;
					time_param[1] = day;
					time_param[2] = hour;
					time_param[3] = min;
					sunray_dir = sunray.changeSunRay(time_param);
					solar_scene->changeSolarScene(sunray_dir);
					double res = sdbk_calc->calcTotalEnergy(sunray.current_DNI);
					outFile << setprecision(12) << sunray.current_altitude << ' ' << sunray.current_azimuth << ' ' << res << endl;
				}
			}
		}
	}
	outFile.close();


	time_param[0] = 1;
	time_param[1] = 1;
	time_param[2] = 12;
	time_param[3] = 0;

	Timer::resetStart();
	sunray_dir = sunray.changeSunRay(time_param);
	solar_scene->changeSolarScene(sunray_dir);
	Timer::printDuration("Change Solar Scene");


	if (options == "-s_l") {

		Timer::resetStart();
		sdbk_calc->calcSampleEnergy(100, 100, sunray.current_DNI);
		Timer::printDuration("Calculate Sample Heliostats' energy");

		Timer::resetStart();
		LSF lsf(50, 50, 100, 100);
		lsf.LSF_surface(solar_scene);
		Timer::printDuration("LSF fitting");
		sdbk_calc->saveCalcRes("lsf_fitting_res.txt");

		Timer::resetStart();
		sdbk_calc->calcTotalEnergy(sunray.current_DNI);
		Timer::printDuration("Calculate total heliostats");
		sdbk_calc->saveCalcRes("total_res.txt");
	}
	

	return 1;
}