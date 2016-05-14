/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include "mirtk/MeanShift.h"

namespace mirtk {

MeanShift::MeanShift(GreyImage& image, int padding, int nBins)
{
	_image = image;
	_orig_image=image;
	_nBins=nBins;
	_padding=padding;
	_brain=NULL;
	_density = new double[_nBins];
	for(int i=0;i<_nBins;i++) _density[i]=0;
	_bg=-1;
	_gm=-1;
	_wm=-1;
	_bin_width=1;

}

MeanShift::~MeanShift()
{
	delete[] _density;
}

RealImage MeanShift::ReturnMask()
{
	RealImage mask = _image;
	RealPixel *ptr = mask.GetPointerToVoxels();
	int n = _image.GetNumberOfVoxels();
	for(int i=0;i<n;i++)
	{
		if(*ptr!=_padding)
			*ptr=1;
		else
			*ptr=0;
		ptr++;
	}
	return mask;
}


void MeanShift::SetOutput(GreyImage *output)
{
	_output = output;
}

double MeanShift::ValueToBin(double value)
{
	return _nBins*(value-_imin)/(_imax+1-_imin);
}

double MeanShift::BinToValue(int bin)
{
	return _imin+bin*(_imax+1-_imin)/_nBins;
}

double MeanShift::GenerateDensity(double cut_off)
{
	int i,j;
	GreyPixel *ptr=_image.GetPointerToVoxels();
	int n = _image.GetNumberOfVoxels();
	int voxels=0;

	_image.GetMinMax(&_imin,&_imax);
	std::cout<<_imin<<" "<<_imax<<std::endl;
	_bin_width = (_imax-_imin+1)/_nBins;

	std::cout<<"generating density..."<<std::endl;
	for(i=0;i<n;i++)
	{
		if(*ptr>_padding)
		{
			j = ValueToBin(*ptr);
			_density[j]++;
			voxels++;
		}
		ptr++;
	}
	std::cout<<"done"<<std::endl;


	std::cout<<"Cutting off 2% of highest intensities"<<std::endl;
	double sum=0,bs=0;
	for(i=0; i<_nBins; i++) sum+=_density[i];
	i=_nBins-1;
	while(bs<cut_off*sum)
	{
		bs+=_density[i];
		i--;
	}

	_limit = BinToValue(i);

	std::cout<<"limit="<<_limit<<std::endl<<std::endl;
	return _limit;
}


void MeanShift::AddPoint(int x, int y, int z)
{
	if((x>=0)&&(y>=0)&&(z>=0)&&(x<_image.GetX())&&(y<_image.GetY())&&(z<_image.GetZ()))
	{
		bool mask=true;
		Point p;
		if(_brain!=NULL) if(_brain->Get(x,y,z)==1) mask=false;

		if((_map.Get(x,y,z)==1)&&(_image.Get(x,y,z)<_treshold)&&(mask))
		{
			p._x=x;
			p._y=y;
			p._z=z;

			_q.push(p);

			_map.Put(x,y,z,0);
		}
	}
}

void MeanShift::AddPoint(int x, int y, int z, int label)
{
	if((x>=0)&&(y>=0)&&(z>=0)&&(x<_image.GetX())&&(y<_image.GetY())&&(z<_image.GetZ()))
	{
		Point p;

		if((_map.Get(x,y,z)==0)&&(_image.Get(x,y,z)==label))
		{
			p._x=x;
			p._y=y;
			p._z=z;

			_q.push(p);

			_map.Put(x,y,z,1);
			_clusterSize++;
		}
	}
}

double MeanShift::msh(double y, double h)
{
	double points, sum, y0;
	GreyPixel* ptr;
	int n = _image.GetNumberOfVoxels();

	//std::cout<<"msh: position = "<<y<<", bandwidth = "<<h<<std::endl;

	if(h<=_bin_width) return y;

	do
	{
		ptr = _image.GetPointerToVoxels();
		sum=0;
		points=0;

		for(int i = 0; i < n; i++)
		{
			if (*ptr > _padding)
			{
				if ( abs(*ptr-y) <= h)
				{
					points++;
					sum += *ptr;
				}
			}
			ptr++;
		}
		y0=y;
		y=sum/points;
	}
	while (abs(y-y0) >= 1);

	return y;
}

double MeanShift::findGMvar()
{
	double k;
	int mean = ValueToBin(_gm);
	double sigma=0;

	for (int i=mean; i>=0; i--)
	{
		k=sqrt(2*log(_density[mean]/(double)_density[i]));
		sigma=abs(BinToValue(i)-_gm)/k;
		std::cout<<"Pos "<<BinToValue(i)<<": sigma = "<<sigma<<std::endl;
	}

	return sigma;
}

double MeanShift::findMax(double tr1, double tr2)
{
	double mx=0, pos=-1;
	int j1,j2,i0,i;
	//j1 = nBins*(tr1-imin)/(imax+1-imin);
	//j2 = nBins*(tr2-imin)/(imax+1-imin);
	j1=ValueToBin(tr1);
	j2=ValueToBin(tr2);

	//std::cout<<imin<<" "<<imax<<std::endl;
	//std::cout<<j1<<" "<<j2<<std::endl;

	i0 = j1;
	for (i=j1; i<=j2; i++)
	{
		//std::cout<<_density[i]<<" ";
		if(mx<_density[i])
		{
			i0=i;
			mx=_density[i];
			//std::cout<<"density["<<i<<"]="<<_density[i]<<", mx="<<mx<<std::endl;
		}
	}
	//std::cout<<i0<<std::endl;
	pos=BinToValue(i0);//imin+i0*(imax+1-imin)/nBins;
	//std::cout<<pos<<std::endl;

	std::cout<<"Maximum for interval <"<<tr1<<", "<<tr2<<"> is "<<pos<<std::endl;
	return pos;
}

double MeanShift::findMin(double tr1, double tr2)
{
	double mx=100000000, pos=-1;
	int j1,j2,i0,i;
	//j1 = nBins*(tr1-imin)/(imax+1-imin);
	//j2 = nBins*(tr2-imin)/(imax+1-imin);

	//std::cout<<imin<<" "<<imax<<std::endl;
	//std::cout<<j1<<" "<<j2<<std::endl;

	j1=ValueToBin(tr1);
	j2=ValueToBin(tr2);
	i0=j1;

	for (i=j1; i<=j2; i++)
	{
		if(mx>_density[i])
		{
			i0=i;
			mx=_density[i];
			//std::cout<<"density["<<i<<"]="<<density[i]<<", mx="<<mx<<std::endl;
		}
	}
	pos=BinToValue(i0);//imin+i0*(imax+1-imin)/nBins;
	//pos=imin+i0*(imax+1-imin)/nBins;

	std::cout<<"Minimum for interval <"<<tr1<<", "<<tr2<<"> is "<<pos<<std::endl;
	return pos;
}

double MeanShift::split(double pos1, double pos2, double bw, double h1, double h2)
{
	double tr=0.02*_nBins;
	double m1, m2, m3;
	bool change = true;

	std::cout<<"split: "<<pos1<<" "<<pos2<<" "<<bw<<" "<<h1<<" "<<h2<<std::endl;

	while(((pos2-pos1) > _bin_width)&&(change)&&( bw>_bin_width ))
	{
		change=false;
		std::cout<<"pos1="<<pos1<<" pos2="<<pos2<<std::endl;
		tr=_bin_width;
		if (tr<1) tr=1;
		m1=msh(pos1,bw);
		m2=msh(pos2,bw);
		std::cout<<"m1="<<m1<<", m2="<<m2<<std::endl;

		if((m2-m1)< tr)
		{
			h2=bw;
			bw=(h1+h2)/2;
			if(bw<=h2-1)
			{
				std::cout<<"reducing bandwidth: bw="<<bw<<std::endl;
				std::cout<<"h1="<<h1<<", h2="<<h2<<std::endl;
				change=true;
			}
		}
		else
		{
			std::cout<<"changing positions:"<<std::endl;
			m3=msh((pos1+pos2)/2, bw);
			std::cout<<"m3="<<m3<<std::endl;
			if((m3-m1)<tr)
			{
				pos1=(pos1+pos2)/2;
				change=true;
				std::cout<<"pos1="<<pos1<<std::endl;
			}
			else
			{
				if((m2-m3)<tr)
				{
					pos2=(pos1+pos2)/2;
					change=true;
					std::cout<<"pos2="<<pos2<<std::endl;
				}
				else
				{
					h1=bw;
					bw=(h1+h2)/2;
					if(bw>=h1+1)
					{
						std::cout<<"increasing bandwidth: bw="<<bw<<std::endl;
						std::cout<<"h1="<<h1<<", h2="<<h2<<std::endl;

						change=true;
					}
					else 
					{
						//std::cout<<"bandwidth fixed:"<<bw<<std::endl;
						std::cout<<"change=false, exiting split."<<std::endl;

					}
				}
			}
			///////end changing positions
		}
	}

	std::cout<<"treshold="<<pos1<<std::endl;

	return pos1;
}


void MeanShift::SetTreshold(double treshold)
{
	std::cout<<"treshold = "<<treshold<<std::endl;

	_treshold = treshold;
}

void MeanShift::SetTreshold()
{

	double pos1=0, pos2=_limit;
	double bw=2*_nBins;
	double h1=0, h2=bw;

	std::cout<<"calculating treshold ... ";

	_limit1 = split(pos1, pos2, bw, h1, h2);
	_bg=findMax(0,_limit1);
	_gm=findMax(_limit1,_limit);
	_treshold = findMin(_bg,_gm);
	//_treshold = _limit1;
}


int MeanShift::Lcc(int label, bool add_second)
{
	int i,j,k;
	int lcc_size = 0;
	int lcc2_size = 0;
	int lcc_x = 0, lcc_y = 0, lcc_z = 0, lcc2_x = 0, lcc2_y = 0, lcc2_z = 0;

	//std::cout<<"Finding Lcc"<<std::endl;
	_map=_image;
	//_image.Write("image.nii.gz");

	GreyPixel* ptr=_map.GetPointerToVoxels();
	int n = _image.GetNumberOfVoxels();
	for(i=0;i<n;i++)
	{
		*ptr=0;
		ptr++;
	}

	for (i=0; i<_image.GetX();i++)
		for (j=0; j<_image.GetY();j++)
			for (k=0; k<_image.GetZ();k++)
			{
				_clusterSize = 0;
				if ((_map.Get(i,j,k)==0)&&(_image.Get(i,j,k)==label)) Grow(i,j,k,label);
				if (_clusterSize > lcc_size)
				{
					lcc2_size = lcc_size;
					lcc2_x = lcc_x;
					lcc2_y = lcc_y;
					lcc2_z = lcc_z;

					lcc_size = _clusterSize;
					lcc_x = i;
					lcc_y = j;
					lcc_z = k;
				}
			}

	ptr=_map.GetPointerToVoxels();
	n = _image.GetNumberOfVoxels();
	for(i=0;i<n;i++)
	{
		*ptr=0;
		ptr++;
	}

	Grow(lcc_x,lcc_y,lcc_z,label);

	if((add_second)&&(lcc2_size > 0.5*lcc_size))
	{
		std::cout<<"Adding second largest cluster too. ";
		Grow(lcc2_x,lcc2_y,lcc2_z,label);
	}
	//_map.Write("lcc.nii.gz");
	*_output = _map;

	return lcc_size;
}

int MeanShift::LccS(int label, double treshold)
{
	int i,j,k;
	int lcc_size = 0;
	int lcc_x, lcc_y, lcc_z;
	queue<Point> seed;
	queue<int> size;

	//std::cout<<"Finding Lcc and all cluster of 70% of the size of Lcc"<<std::endl;
	_map=_image;
	//_image.Write("image.nii.gz");

	GreyPixel* ptr=_map.GetPointerToVoxels();
	int n = _image.GetNumberOfVoxels();
	for(i=0;i<n;i++)
	{
		*ptr=0;
		ptr++;
	}

	for (i=0; i<_image.GetX();i++)
		for (j=0; j<_image.GetY();j++)
			for (k=0; k<_image.GetZ();k++)
			{
				_clusterSize = 0;
				if ((_map.Get(i,j,k)==0)&&(_image.Get(i,j,k)==label)) Grow(i,j,k,label);
				if (_clusterSize > lcc_size)
				{
					lcc_size = _clusterSize;
					lcc_x = i;
					lcc_y = j;
					lcc_z = k;
				}

				if (_clusterSize > 0)
				{
					size.push(_clusterSize);
					Point p(i,j,k);
					seed.push(p);
				}
			}

	ptr=_map.GetPointerToVoxels();
	n = _image.GetNumberOfVoxels();
	for(i=0;i<n;i++)
	{
		*ptr=0;
		ptr++;
	}

	while (!size.empty())
	{
		if ( size.front()> treshold*lcc_size ) Grow(seed.front()._x,seed.front()._y,seed.front()._z,label);
		seed.pop();
		size.pop();
	}
	//_map.Write("lcc.nii.gz");
	*_output = _map;

	return lcc_size;
}


void MeanShift::Grow(int x, int y, int z, int label)
{
	AddPoint(x,y,z,label);

	while (!_q.empty())
	{
		x=_q.front()._x;
		y=_q.front()._y;
		z=_q.front()._z;

		_q.pop();

		AddPoint(x-1,y,z,label);
		AddPoint(x,y-1,z,label);
		AddPoint(x,y,z-1,label);
		AddPoint(x+1,y,z,label);
		AddPoint(x,y+1,z,label);
		AddPoint(x,y,z+1,label);
	}

}


void MeanShift::RegionGrowing()
{
	int i; 
	std::cout<<"Removing background"<<std::endl;
	_map=_image;

	GreyPixel* ptr=_map.GetPointerToVoxels();
	int n = _image.GetNumberOfVoxels();

	for(i=0;i<n;i++)
	{
		*ptr=1;
		ptr++;
	}

	AddPoint(0,0,0);
	AddPoint(_image.GetX()-1,0,0);
	AddPoint(0,_image.GetY()-1,0);
	AddPoint(0,0,_image.GetZ()-1);
	AddPoint(_image.GetX()-1,_image.GetY()-1,0);
	AddPoint(_image.GetX()-1,0,_image.GetZ()-1);
	AddPoint(0,_image.GetY()-1,_image.GetZ()-1);
	AddPoint(_image.GetX()-1,_image.GetY()-1,_image.GetZ()-1);

	int x,y,z;

	while (!_q.empty())
	{
		x=_q.front()._x;
		y=_q.front()._y;
		z=_q.front()._z;

		_q.pop();

		AddPoint(x-1,y,z);
		AddPoint(x,y-1,z);
		AddPoint(x,y,z-1);
		AddPoint(x+1,y,z);
		AddPoint(x,y+1,z);
		AddPoint(x,y,z+1);
	}
}


void MeanShift::RemoveBackground()
{
	int i; 
	std::cout<<"Removing background"<<std::endl;
	_map=_image;

	GreyPixel* ptr=_map.GetPointerToVoxels();
	int n = _image.GetNumberOfVoxels();

	for(i=0;i<n;i++)
	{
		*ptr=1;
		ptr++;
	}

	AddPoint(0,0,0);
	AddPoint(_image.GetX()-1,0,0);
	AddPoint(0,_image.GetY()-1,0);
	AddPoint(0,0,_image.GetZ()-1);
	AddPoint(_image.GetX()-1,_image.GetY()-1,0);
	AddPoint(_image.GetX()-1,0,_image.GetZ()-1);
	AddPoint(0,_image.GetY()-1,_image.GetZ()-1);
	AddPoint(_image.GetX()-1,_image.GetY()-1,_image.GetZ()-1);

	int x,y,z;

	while (!_q.empty())
	{
		x=_q.front()._x;
		y=_q.front()._y;
		z=_q.front()._z;

		_q.pop();

		AddPoint(x-1,y,z);
		AddPoint(x,y-1,z);
		AddPoint(x,y,z-1);
		AddPoint(x+1,y,z);
		AddPoint(x,y+1,z);
		AddPoint(x,y,z+1);
	}

	std::cout<< "dilating and eroding ... ";
	_brain = new GreyImage(_map);

	int iterations = 3;
	Dilation<GreyPixel> dilation;
	dilation.Input(_brain);
	dilation.Output(_brain);
	for (i = 0; i < iterations; i++) dilation.Run();

	Erosion<GreyPixel> erosion;
	erosion.Input(_brain);
	erosion.Output(_brain);
	for (i = 0; i < iterations; i++) erosion.Run();

	std::cout<<"recalculating ... ";

	ptr=_map.GetPointerToVoxels();
	for(i=0;i<n;i++)
	{
		*ptr=1;
		ptr++;
	}

	AddPoint(0,0,0);
	AddPoint(_image.GetX()-1,0,0);
	AddPoint(0,_image.GetY()-1,0);
	AddPoint(0,0,_image.GetZ()-1);
	AddPoint(_image.GetX()-1,_image.GetY()-1,0);
	AddPoint(_image.GetX()-1,0,_image.GetZ()-1);
	AddPoint(0,_image.GetY()-1,_image.GetZ()-1);
	AddPoint(_image.GetX()-1,_image.GetY()-1,_image.GetZ()-1);

	while (!_q.empty())
	{
		x=_q.front()._x;
		y=_q.front()._y;
		z=_q.front()._z;

		_q.pop();

		AddPoint(x-1,y,z);
		AddPoint(x,y-1,z);
		AddPoint(x,y,z-1);
		AddPoint(x+1,y,z);
		AddPoint(x,y+1,z);
		AddPoint(x,y,z+1);
	}

	std::cout<<"eroding ... ";

	*_brain = _map;

	erosion.Input(_brain);
	erosion.Output(_brain);
	for (i = 0; i < iterations; i++) erosion.Run();

	std::cout<<"final recalculation ...";

	ptr=_map.GetPointerToVoxels();
	for(i=0;i<n;i++)
	{
		*ptr=1;
		ptr++;
	}

	AddPoint(0,0,0);
	AddPoint(_image.GetX()-1,0,0);
	AddPoint(0,_image.GetY()-1,0);
	AddPoint(0,0,_image.GetZ()-1);
	AddPoint(_image.GetX()-1,_image.GetY()-1,0);
	AddPoint(_image.GetX()-1,0,_image.GetZ()-1);
	AddPoint(0,_image.GetY()-1,_image.GetZ()-1);
	AddPoint(_image.GetX()-1,_image.GetY()-1,_image.GetZ()-1);

	while (!_q.empty())
	{
		x=_q.front()._x;
		y=_q.front()._y;
		z=_q.front()._z;

		_q.pop();

		AddPoint(x-1,y,z);
		AddPoint(x,y-1,z);
		AddPoint(x,y,z-1);
		AddPoint(x+1,y,z);
		AddPoint(x,y+1,z);
		AddPoint(x,y,z+1);
	}

	delete _brain;
	std::cout<<"done."<<std::endl;

	ptr=_map.GetPointerToVoxels();
	GreyPixel *ptr_b=_image.GetPointerToVoxels();
	for(i=0;i<n;i++)
	{
		if(*ptr==0) *ptr_b = _padding;
		ptr++;
		ptr_b++;
	}


}

void MeanShift::Write(char *output_name)
{
	_image.Write(output_name);
}

void MeanShift::WriteMap(char *output_name)
{
	_map.Write(output_name);
}

void MeanShift::FindWMGMmeans()
{
	std::cout<<"Finding WM and GM mean"<<std::endl;
	std::cout<<"Make sure that background has been removed when calling this method!!!"<<std::endl;

	_gm=findMax(0, _limit);
	_limit2 = split(_gm, _limit, 2*_nBins, 0, 2*_nBins);
	_wm=findMax(_limit2,_limit);
	_split2=findMin(_gm,_wm);
	std::cout<<"GM mean = "<<_gm<<std::endl;
	std::cout<<"WM mean = "<<_wm<<std::endl;

	/*
  _limit2 = split(0, _limit, 2*_nBins, 0, 2*_nBins);


  _gm=findMax(0, _limit2);
  _wm=findMax(_limit2,_limit);

  if(_bg>-1) _split1=findMin(_bg,_gm);
  _split2=findMin(_gm,_wm);

  std::cout<<"GM mean = "<<_gm<<std::endl;
  std::cout<<"WM mean = "<<_wm<<std::endl;

  GreyImage gmwm(_image);
  GreyPixel *ptr=gmwm.GetPointerToVoxels();
  int n=gmwm.GetNumberOfVoxels();
  for(int i=0;i<n;i++)
  {
    if(*ptr>_padding)
    {
      if(*ptr < _split2) *ptr=0;
      else *ptr=1;
    }
    ptr++;
  }
  gmwm.Write("gm-wm-split.nii.gz");

	 */
}

}



