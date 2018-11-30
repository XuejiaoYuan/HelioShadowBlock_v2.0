//
// Created by Amber on 2018/4/3.
//

#include "SolarScene.h"

//
//  Load heliostat field's file with Duan's format
//
bool SolarScene::initSolarScene(const string &scene_filepath, const Vector3f &sunray_dir){
	this->sunray_dir = Vector3f(sunray_dir.x(), sunray_dir.y(), sunray_dir.z());
	fstream inFile(scene_filepath);
    if(inFile.fail()){
        cerr<<"Can't open the file!"<<endl;
        return false;
    }

    InputMode input_mode = Initial;
    string line, word;
    stringstream line_stream;
    int grid_num;
    int helio_type;
    Receiver* recv;
    Layout* layout;
    Heliostat* helio;
	Vector2f helio_gap;
	Vector2i helio_matrix;
    while(getline(inFile,line)){
        line_stream.clear();
        line_stream.str(line);
        line_stream >> word;

        if(word=="#"){
            line_stream >> word;
            if(word == "Ground"){
                input_mode = GroundMode;
                continue;
            }
            else if(word == "Receiver"){
                input_mode = ReceiverMode;
                continue;
            }
            else if(word == "Heliostats"){
                input_mode = HeliostatMode;
                continue;
            }
            else{
                input_mode = LayoutMode;
                continue;
            }
        }

        switch(input_mode){
            case GroundMode:{
                if(word == "ground")
                    line_stream >> scene_length >>scene_width;
                else if(word == "ngrid"){
                    line_stream >> grid_num;
                    input_mode = Initial;
                }
                break;
            }
            case ReceiverMode:{
                int recv_type;
                ReceiverCreator recv_creator;
                line_stream >> recv_type;
                recv = recv_creator.getReceiver((ReceiverType)recv_type);
				recv->init_recv(inFile, input_mode);
                recvs.push_back(recv);
                recv = nullptr;
                break;
            }
			case HeliostatMode: {
				if (word == "gap")
					line_stream >> helio_gap.x() >> helio_gap.y();
				else if (word == "matrix")
					line_stream >> helio_matrix.x() >> helio_matrix.y();
				else if (word == "end") {
					input_mode = Initial;
					break;
				}
				else if (word == "helio") {
					helio = new Heliostat((HelioType)helio_type);
					helio->initHeliostat(line_stream, inFile, layouts[0]->layout_type, helio_gap, helio_matrix, recvs[0]->focus_center, sunray_dir);
					helio->helio_index = helios.size();
					helios.push_back(helio);
					helio = nullptr;		
				}
				break;
			}
            case LayoutMode:{
                int layout_type;
                LayoutCreator layout_creator;
                line_stream >> layout_type;
                layout = layout_creator.getLayout((LayoutType)layout_type);
				layout->initLayout(inFile, input_mode, helio_type);
                layouts.push_back(layout);
                layout = nullptr;
                break;
            }
            case Initial:break;
            default: break;
        }

    }
    inFile.close();
	layouts[0]->setHelioLayout(helios);
    return true;
}


bool SolarScene::changeSolarScene(const Vector3f & sunray_dir)
{
	this->sunray_dir = Vector3f(sunray_dir.x(), sunray_dir.y(), sunray_dir.z() );

#pragma omp parallel for
	for (int i = 0; i < helios.size(); i++)
		helios[i]->changeSurfaceNormal(recvs[0]->focus_center, sunray_dir);

	layouts[0]->setHelioLayout(helios);
	
	return true;
}

//
// [镜场优化预处理] 镜场优化中固定参数设置
//
bool SolarScene::initFieldParam(const string& file_name)
{
	fstream inFile(file_name, ios_base::in);
	if (inFile.fail()) {
		cerr << "Can't open the filed parameter file!" << endl;
		return false;
	}

	string line, word;
	stringstream line_stream;
	InputMode input_mode = Initial;
	int row = 0;
	int col = 0;
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
			recvs.push_back(recv);
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
			else if (word == "row")
				line_stream >> row;
			else if (word == "col")
				line_stream >> col;
			else if (word == "type")
				line_stream >> helio_type;
			else if (word == "end") {
				input_mode = Initial;
				layouts.push_back(layout);
				layout = nullptr;
			}
			break;
		}
		case Initial:
			break;
		default:
			break;
		}
	}
	inFile.close();
	layouts[0]->helio_type = (HelioType)helio_type;
	layouts[0]->helio_pos = helio_pos;
	layouts[0]->helio_gap = helio_gap;
	layouts[0]->helio_matrix = helio_matrix;
	layouts[0]->helio_size = helio_size;

	return true;
}

//
// [镜场优化] 镜场优化时改变优化参数以调整镜场位置
//
bool SolarScene::adjustFieldParam(const vector<vector<float>*>& field_args)
{
	if (!helios.empty()) {
		for (auto&h : helios)
			delete h;
	}
	helios.clear();
	layouts[0]->adjustHelioLayout(helios, field_args, recvs);
	return true;
}
