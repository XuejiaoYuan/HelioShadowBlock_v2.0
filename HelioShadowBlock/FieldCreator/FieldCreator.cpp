#include "FieldCreator.h"

void FieldCreator::initFiledParam(string & file_name, SolarScene * scene)
{
	fstream inFile(file_name, ios_base::in);
	if (inFile.fail()) {
		cerr << "Can't open the filed parameter file!" << endl;
	}
	if (scene != NULL)
		delete scene;
	scene = new SolarScene;

	string line, word;
	stringstream line_stream;
	InputMode input_mode = Initial;
	Heliostat* helio;
	int nhelio = 0;
	Layout* layout;
	int helio_type;
	int layout_type;
	Vector2f helio_gap;
	Vector2i helio_matrix;
	Vector3f helio_size;
	Vector3f helio_pos;
	while (getline(inFile, line)) {
		line_stream.clear();
		line_stream.str(line);
		line_stream >> word;

		if (word == "#") {
			line_stream >> word;
			if (word == "Receiver") {
				input_mode = ReceiverMode;
				continue;
			}
			else if (word == "Heliostats") {
				input_mode = HeliostatMode;
				continue;
			}
			else if (word == "Grid") {
				input_mode = LayoutMode;
				continue;
			}
		}

		switch (input_mode)
		{
		case ReceiverMode: {
			int recv_type;
			ReceiverCreator recv_creator;
			line_stream >> recv_type;
			Receiver* recv = recv_creator.getReceiver((ReceiverType)recv_type);
			recv->init_recv(inFile, input_mode);
			scene->recvs.push_back(recv);
			break;
		}
		case HeliostatMode: {
			if (word == "gap")
				line_stream >> helio_gap.x() >> helio_gap.y();
			else if (word == "matrix")
				line_stream >> helio_matrix.x() >> helio_matrix.y();
			else if (word == "pos")
				line_stream >> helio_pos.x() >> helio_pos.y() >> helio_pos.z();
			else if (word == "size")
				line_stream >> helio_size.x() >> helio_size.y() >> helio_size.z();
			else if (word == "end")
				input_mode = Initial;
			break;
		}
		case LayoutMode: {
			if (word == "Scene") {
				line_stream >> layout_type;
				LayoutCreator layout_creator;
				layout = layout_creator.getLayout((LayoutType)layout_type);
			}
			else if (word == "nhelio")
				line_stream >> nhelio;
			else if (word == "type")
				line_stream >> helio_type;
			else if (word == "end")
				input_mode = Initial;
			scene->layouts.push_back(layout);
			layout = nullptr;
			break;
		}
		case Initial:
			break;
		default:
			break;
		}
	}
	inFile.close();
	for (int i = 0; i < nhelio; i++) {
		helio = new Heliostat((HelioType)helio_type);
		helio->helio_gap = helio_gap;
		helio->helio_matrix = helio_matrix;
		helio->helio_size = helio_size;
		helio->helio_pos = helio_pos;
		scene->helios.push_back(helio);
		helio = nullptr;
	}

	//layouts[0]->setHelioLayout(helios);
	return;
}
