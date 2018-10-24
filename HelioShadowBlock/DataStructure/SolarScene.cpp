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
    //HeliostatCreator helio_creator;
	Vector3f focus_center = Vector3f(0, 0, 0);
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
                while(getline(inFile, line)) {
                    line_stream.clear();
                    line_stream.str(line);
                    line_stream >> word;
                    if (word == "end"){
                        input_mode = Initial;
                        break;
                    }
                    else if (word == "pos")
                        line_stream >> recv->recv_pos.x() >> recv->recv_pos.y() >> recv->recv_pos.z();
                    else if (word == "size")
                        line_stream >> recv->recv_size.x() >> recv->recv_size.y() >> recv->recv_size.z();
                    else if (word == "norm")
                        line_stream >> recv->recv_normal.x() >> recv->recv_normal.y() >> recv->recv_normal.z();
                    else
                        line_stream >> recv->recv_face;

                }
                focus_center = recv->recv_pos + Vector3f(recv->recv_normal.array() * recv->recv_size.array());
                recvs.push_back(recv);
                recv = nullptr;
                break;
            }
            case HeliostatMode:{
                if(word == "gap")
                    line_stream >> helio_gap.x() >> helio_gap.y();
                else if(word == "matrix")
                    line_stream >> helio_matrix.x() >> helio_matrix.y();
                else if(word == "helio"){
                    //helio = helio_creator.getHeliostat((HelioType)helio_type);
					helio = new Heliostat((HelioType)helio_type);
					line_stream >> helio->helio_pos.x() >> helio->helio_pos.y() >> helio->helio_pos.z() ;
					getline(inFile, line);
					line_stream.clear();
					line_stream.str(line);
                    line_stream >> helio->helio_size.x() >> helio->helio_size.y() >> helio->helio_size.z() ;
                    helio->helio_gap = helio_gap;
                    helio->helio_matrix = helio_matrix;
                    bool flag = helio->initSurfaceNormal(focus_center, sunray_dir);
                    if(!flag)
                        return false;
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
	// layouts[0]->setHelioLayout(helios);
	for (auto& layout : layouts)
		layout->setHelioLayout(helios);
    return true;
}

bool SolarScene::changeSolarScene(const Vector3f & sunray_dir)
{
	this->sunray_dir = Vector3f(sunray_dir.x(), sunray_dir.y(), sunray_dir.z() );
	Vector3f focus_center = recvs[0]->recv_pos + Vector3f(recvs[0]->recv_normal.array() * recvs[0]->recv_size.array());

	for (auto&helio : helios) {
		helio->changeSurfaceNormal(focus_center, sunray_dir);
	}
	return true;
}

