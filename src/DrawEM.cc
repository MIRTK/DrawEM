/*
 * Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Christian Ledig
 * Copyright 2013-2016 Antonios Makropoulos
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


#include "mirtk/DrawEM.h"

namespace mirtk {

// Default constructor
DrawEM::DrawEM() : EMBase(){
    InitialiseParameters();
}

template <class ImageType>
DrawEM::DrawEM(int noTissues, ImageType **atlas, ImageType *background) : EMBase(noTissues, atlas, background)
{
    InitialiseParameters();
}

template <class ImageType>
DrawEM::DrawEM(int noTissues, ImageType **atlas ) : EMBase(noTissues, atlas )
{
    InitialiseParameters();
}

template <class ImageType>
DrawEM::DrawEM(int noTissues, ImageType **atlas, ImageType **initposteriors ) : EMBase(noTissues, atlas, initposteriors )
{
    InitialiseParameters();
}

void DrawEM::InitialiseParameters(){
    outlabel=1;
    csflabel=2;
    gmlabel=3;
    wmlabel=4;
    mrfweight=1;
    bignn=false;
}


void DrawEM::SetInput(const RealImage &image, const Matrix &connectivity)
{
    _uncorrected = image;
    if( connectivity.Cols() != _number_of_tissues || connectivity.Rows() != connectivity.Cols() )
    {
        std::cerr << "Warning: Connectivity matrix has wrong size! expected: " << _number_of_tissues << "x"<< _number_of_tissues << std::endl;
        std::cerr << "is:" << connectivity.Rows() << "x" << connectivity.Rows() << std::endl;
    }
    _connectivity = connectivity;

    EMBase::SetInput(_uncorrected);


    //SEG {
    huipvcorr=false;
    intermrf=false;
    beta = 1.0/3.0;
    betainter = 1.0/3.0;
    //SEG }
}


void DrawEM::BStep()
{
    // Create bias correction filter
    _biascorrection.SetInput(&_uncorrected, &_estimate);
    _biascorrection.SetWeights(&_weights);
    _biascorrection.SetOutput(_biasfield);
    _biascorrection.SetPadding((short int) _padding);
    _biascorrection.SetMask(&_mask);
    _biascorrection.Run();

    // Generate bias corrected image for next iteration
    _input = _uncorrected;
    _biascorrection.Apply(_input);
}


// Relaxation according to Cardoso in MICCAI 2011
void DrawEM::RStep(){
    RStep(0.5);
}

void DrawEM::RStep(double rf)
{
    double relaxFactor = rf;

    _atlas.First();
    _output.First();
    HashProbabilisticAtlas filteredAtlas;

    RealImage filterInput;
    for( int k = 0; k < _number_of_tissues; ++k )
    {
        filterInput = _output.GetImage(k).ToGenericImage();
        GaussianBlurring<RealPixel> filter(2.0);
        filter.Input(&filterInput);
        filter.Output(&filterInput);
        filter.Run();
        filteredAtlas.AddImage(filterInput);
    }

    filteredAtlas.First();
    RealPixel *ptr = _input.GetPointerToVoxels();
    BytePixel *pm = _mask.GetPointerToVoxels();


    //mine
    bool bMRF = _number_of_tissues == _connectivity.Rows();



    int per = 0;
    double* numerator = new double[_number_of_tissues];
    double values[_number_of_tissues];

    for (int i=0; i< _number_of_voxels; i++){
        if (i*10.0/_number_of_voxels > per) {
            per++;
            std::cerr<<per<<"0%...";
        }
        double denominator = 0.0;
        if (*pm == 1){
            for( int k = 0; k < _number_of_tissues; ++k ){
                values[k]= (1.0-relaxFactor) * filteredAtlas.GetValue(k) + relaxFactor * _atlas.GetValue(k);
                filteredAtlas.SetValue(k, values[k]);
            }

            for( int k = 0; k < _number_of_tissues; ++k ){
                numerator[k] = .0;
                double temp = .0;

                if(bMRF){
                    for( int j = 0; j < _number_of_tissues; ++j ){
                        if ( _connectivity(k,j) <= 1 ){
                            temp += values[k] * values[j];
                        }
                    }
                }else{
                    temp = values[k];
                }
                numerator[k] += temp;
                denominator += numerator[k];
            }

            for( int k = 0; k < _number_of_tissues; ++k ){
                if( denominator != 0 ){
                    _atlas.SetValue(k, numerator[k] / denominator);
                }else{
                    _atlas.SetValue(k, 1/k);
                }
            }

            if (denominator <= 0) {
                int x,y,z;
                _input.IndexToVoxel(i, x, y, z);
                std::cerr<<"Division by 0 while computing relaxed prior probabilities at voxel "<<x<<","<<y<<","<<z<<std::endl;
            }
        }
        pm++;
        ptr++;
        _output.Next();
        _atlas.Next();
        filteredAtlas.Next();
    }
    delete[] numerator;
}





int DrawEM::AddPartialVolumeClass(int classA, int classB, int huiclass)
{
    RealImage pvclass = _input;
    RealPixel *ptr_pvclass = pvclass.GetPointerToVoxels();

    // mixing coefficient
    double gamma = 0.0;
    int N = 0;
    RealPixel *ptr = _input.GetPointerToVoxels();
    BytePixel *pm = _mask.GetPointerToVoxels();

    for( int i = 0; i < pvclass.GetNumberOfVoxels(); ++i )
    {
        if( *pm == 1 )
        {
            double fc = (_mi[classA] - *ptr ) / ( _mi[classA]-_mi[classB]);
            if( fc >= 0 && fc <= 1.0 )
            {
                gamma += fc;
                N++;
            }
        }
        *ptr_pvclass = .0;
        ptr_pvclass++;
        pm++;
        ptr++;
    }
    if( N > 0 )
    {
        gamma /= N;
    }
    else
    {
        std::cerr << "No mixel voxels found, not adding partial volume class!" << std::endl;
        return -1;
    }

    double* mi = new double[_number_of_tissues+1];
    double* sigma = new double[_number_of_tissues+1];

    mi[_number_of_tissues] = 0;
    sigma[_number_of_tissues] = 0;
    for( int i = 0; i < _number_of_tissues; ++i )
    {
        mi[i] = _mi[i];
        sigma[i] = _sigma[i];
    }
    mi[_number_of_tissues] = (1.0 - gamma) * mi[classA] + gamma * mi[classB];
    sigma[_number_of_tissues] = (1.0 - gamma) * (1.0 - gamma) * sigma[classA] + gamma * gamma * sigma[classB];

    _atlas.First();
    _output.First();
    ptr = pvclass.GetPointerToVoxels();
    pm = _mask.GetPointerToVoxels();

    for( int i = 0; i < _number_of_voxels; ++i )
    {
        if( *pm == 1 )
        {
            double tmp = _output.GetValue(classA) * _output.GetValue(classB);
            if( tmp > 0.0 )
            {
                //add new PV class with ω∗i(j/k) =√pij pik

                tmp = sqrt(tmp) / 0.5;
            }
            else
            {
                tmp = 0.0;
            }
            *ptr = tmp;
            for( int k = 0; k < _number_of_tissues; ++k )
            {
                tmp += _output.GetValue(k);
            }
            if( tmp > 0.0 )
            {
                *ptr /= tmp;
                for( int k = 0; k < _number_of_tissues; ++k )
                {
                    _atlas.SetValue(k,_output.GetValue(k) / tmp);
                }
            }
            else
            {
                std::cerr << "error probability = 0" << std::endl;
                return -1;
            }

        }
        ptr++;
        pm++;
        _atlas.Next();
        _output.Next();
    }

    delete[] _mi;
    delete[] _sigma;

    _number_of_tissues++;
    _mi = new double[_number_of_tissues];
    _sigma = new double[_number_of_tissues];
    for( int k = 0; k < _number_of_tissues; ++k)
    {
        _mi[k] = mi[k];
        _sigma[k] = sigma[k];
    }
    delete[] mi;
    delete[] sigma;

    _atlas.AddImage(pvclass);

    RealImage newimage = pvclass;
    _output.AddImage(newimage);

    _atlas.First();
    _output.First();
    for( int i = 0; i < _number_of_voxels; ++i )
    {
        for( int k = 0; k < _number_of_tissues; ++k)
        {
            _output.SetValue(k, _atlas.GetValue(k));
        }
        _atlas.Next();
        _output.Next();
    }
    std::cout << "Connectivity before update" << std::endl;
    _connectivity.Print();
    // PV classes get same connectivity as parent class!
    Matrix newconnectivity(_connectivity.Rows()+1, _connectivity.Cols()+1);
    for( int i = 0; i < newconnectivity.Rows(); ++i )
    {
        for( int j = 0; j < newconnectivity.Cols(); ++j )
        {
            if( i < _connectivity.Rows() && j < _connectivity.Cols() )
            {
                // classA and classB are now distant to each other!
                if( (i == classA && j == classB) || (i == classB && j == classA) )
                {
                    newconnectivity.Put(i,j,2);
                }
                // rest stays as it was
                else
                {
                    newconnectivity.Put(i,j,_connectivity.Get(i,j));
                }
            }
            // same class
            else if( i == j )
            {
                newconnectivity.Put(i,j,0);
            }
            // this is for the PV class, close to pv contributing classes, distant to all others
            else
            {
                if( i == classA || i == classB || j == classA || j == classB )
                {
                    newconnectivity.Put(i,j, 1);
                }
                else
                {
                    newconnectivity.Put(i,j, 2);
                }
            }
        }
    }

    // if classification is done with backgound ensure that background class is the last tissue class
    // -> swap prob maps and adapt connectivity matrix
    if( _has_background )
    {
        double tmp =_mi[_number_of_tissues-1];
        _mi[_number_of_tissues-1] = _mi[_number_of_tissues-2];
        _mi[_number_of_tissues-2] = tmp;
        _sigma[_number_of_tissues-1] = _sigma[_number_of_tissues-2];
        _sigma[_number_of_tissues-2] = tmp;

        //_output.SwapImages(_number_of_tissues-1, _number_of_tissues-2);
        //_atlas.SwapImages(_number_of_tissues-1, _number_of_tissues-2);

        // now swap last two rows and columns of newconnectivity, since prob maps (for background) were swapped as well
        for( int i = 0; i < newconnectivity.Cols(); ++i )
        {
            double tmp = newconnectivity.Get(i,newconnectivity.Rows()-1);
            newconnectivity.Put(i,newconnectivity.Rows()-1, newconnectivity.Get(i,newconnectivity.Rows()-2));
            newconnectivity.Put(i,newconnectivity.Rows()-2, tmp);
        }
        for( int i = 0; i < newconnectivity.Rows(); ++i )
        {
            double tmp = newconnectivity.Get(newconnectivity.Cols()-1,i);
            newconnectivity.Put(newconnectivity.Cols()-1,i, newconnectivity.Get(newconnectivity.Cols()-2,i));
            newconnectivity.Put(newconnectivity.Cols()-2, i, tmp);
        }
    }

    _connectivity = newconnectivity;
    int pv_position = _number_of_tissues - 1;
    //mine
    //int pv_position = _number_of_tissues;
    if( _has_background )
    {
        pv_position = _number_of_tissues - 2;
        //mine
        //pv_position = _number_of_tissues - 1;
    }
    pv_classes.insert(make_pair(pv_position, pv_connections.size() ) );

    pv_connections.push_back( make_pair(classA, classB) );
    pv_fc.push_back(gamma);
    std::cout << "Fractional Content of classA=" << 1.0-gamma << std::endl;
    std::cout << "connectivity after update " << std::endl;
    _connectivity.Print();


    int *pretissuelabels=new int[_number_of_tissues-1];
    for(int i=0;i<_number_of_tissues-1;i++){
        pretissuelabels[i]=tissuelabels[i];
    }

    delete tissuelabels;
    tissuelabels=new int[_number_of_tissues];
    for(int i=0;i<_number_of_tissues-1;i++){
        tissuelabels[i]=pretissuelabels[i];
    }
    tissuelabels[_number_of_tissues-1]=huiclass;



    return pv_position;
}


bool print=true;
double DrawEM::getMRFenergy(int index, int tissue)
{


    if( _connectivity.Rows() == 1 )
    {
        return 1.0;
    }
    double dx,dy,dz;
    _input.GetPixelSize(&dx,&dy,&dz);

    double sx, sy,sz;
    sx = 1.0/dx;
    sy = 1.0/dy;
    sz = 1.0/dz;

    double expo = -1.0* mrfweight;

    double energy = .0;
    int x,y,z;
    _input.IndexToVoxel(index, x, y, z);
    expo *= beta;

    double weight = 0;
    int lx = max(x-1,0);
    int rx = min(x+1, _input.GetX()-1 );
    int ly = max(y-1,0);
    int ry = min(y+1, _input.GetY()-1 );
    int lz = max(z-1,0);
    int rz = min(z+1, _input.GetZ()-1 );

    for( int k = 0; k < _number_of_tissues; k++)
    {
        double temp = 0;
        weight = 0;
        double conn = _connectivity.Get(k, tissue);
        weight=conn;
        if( conn == 2 ) weight = 5;
        else if (conn==0) continue;


        temp += sx * (_output.GetValue( rx ,y,z,k) + _output.GetValue( lx,y,z,k) );
        temp += sy * ( _output.GetValue(x,ry,z,k) + _output.GetValue(x,ly,z,k) );
        temp += sz * ( _output.GetValue(x,y,rz,k) + _output.GetValue(x,y,lz,k) );

        energy += temp * weight;
    }

    expo *= energy;
    return exp(expo);
}



double wneighbors[3][3][3];
double DrawEM::getMRFenergy_diag(int index, int tissue)
{


    if( _connectivity.Rows() == 1 )
    {
        return 1.0;
    }
    double dx,dy,dz;
    _input.GetPixelSize(&dx,&dy,&dz);



    double expo = -1.0;

    double energy = .0;
    int x,y,z;
    _input.IndexToVoxel(index, x, y, z);
    expo *= beta;

    double weight = 0;

    int lx = max(x-1,0);
    int rx = min(x+1, _input.GetX()-1 );
    int ly = max(y-1,0);
    int ry = min(y+1, _input.GetY()-1 );
    int lz = max(z-1,0);
    int rz = min(z+1, _input.GetZ()-1 );

    double dist;
    for(int cx=lx;cx<=rx;cx++){
        for(int cy=ly;cy<=ry;cy++){
            for(int cz=lz;cz<=rz;cz++){
                if(cx==x && cy==y && cz==z) continue;
                dist=1.0/sqrt(pow(dx*(cx-x),2)+pow(dy*(cy-y),2)+pow(dz*(cz-z),2) );
                wneighbors[cx-x+1][cy-y+1][cz-z+1]=dist;
            }
        }
    }



    for( int k = 0; k < _number_of_tissues; k++)
    {
        double temp = 0;
        weight = 0;
        double conn = _connectivity.Get(k, tissue);

        weight=conn;

        for(int cx=lx;cx<=rx;cx++)
            for(int cy=ly;cy<=ry;cy++)
                for(int cz=lz;cz<=rz;cz++){
                    if(cx==x && cy==y && cz==z) continue;
                    temp += wneighbors[cx-x+1][cy-y+1][cz-z+1]*_output.GetValue( cx ,cy,cz,k);
                }

        energy += temp * weight;//    	G.Get(k,tissue);
    }

    expo *= energy;
    expo /=2;

    return exp(expo);
}




void DrawEM::EStepMRF()
{
    std::cout << "E-step with MRF" <<std::endl;

    IntegerImage segmentation;

    int i, k;
    double x;
    Gaussian* G = new Gaussian[_number_of_tissues];

    for (k = 0; k < _number_of_tissues; k++) {
        G[k].Initialise( _mi[k], _sigma[k]);
    }

    RealPixel *pptr;
    if(_postpen) pptr= _postpenalty.GetPointerToVoxels();

    _atlas.First();
    _output.First();
    RealPixel *ptr = _input.GetPointerToVoxels();
    BytePixel *pm = _mask.GetPointerToVoxels();
    int per = 0;
    double* numerator = new double[_number_of_tissues];
    double denominator=0, temp=0;
    bool bMRF = _number_of_tissues == _connectivity.Rows();


    for (i=0; i< _number_of_voxels; i++) {
        if (i*10.0/_number_of_voxels > per) {
            per++;
            std::cout<<per<<"0%...";
        }
        denominator = 0;
        temp = 0;

        if (*pm == 1) {


            x = *ptr;
            double MRFenergies[_number_of_tissues];
            double denominatorMRF = .0;

            for (k = 0; k < _number_of_tissues; k++) {
                double mrfenergy=1;

                if(beta!=0){
                    if(bignn)mrfenergy*=getMRFenergy_diag(i,k);
                    else mrfenergy*=getMRFenergy(i,k);
                }

                if(intermrf && betainter!=0) mrfenergy*=getMRFInterEnergy(i,k);
                MRFenergies[k] = _atlas.GetValue(k) * mrfenergy;
                denominatorMRF += MRFenergies[k];
            }




            for (k = 0; k < _number_of_tissues; k++) {
                temp = G[k].Evaluate(x);

                // MRF matrix fits number of tissues?
                if( bMRF )
                {
                    //different step of MRF!!
                    temp = temp * MRFenergies[k] / denominatorMRF;
                }
                else
                {
                    temp = temp * _atlas.GetValue(k);
                }

                numerator[k] = temp;
                denominator += temp;
            }

            //model averaging
            if (_postpen && denominator != 0) {
                double olddenom=denominator;
                denominator=0;
                for (k = 0; k < _number_of_tissues; k++) {
                    double value = numerator[k]/olddenom;
                    double priorvalue=_atlas.GetValue(k);
                    value=(1-*pptr)*value +*pptr *  priorvalue;
                    numerator[k] = value;
                    denominator += value;
                }
            }



            for (k = 0; k < _number_of_tissues; k++) {
                if (denominator > 0) {
                    double value = numerator[k]/denominator;
                    if ((value < 0) || (value > 1)) {
                        int x,y,z;
                        _input.IndexToVoxel(i, x, y, z);
                        std::cerr << "Probability value = " << value <<" @ Estep-mrf at voxel "<< x<<" "<<y<<" "<<z<< ", structure " << k << std::endl;
                        //exit(1);
                        if (value < 0)value=0;
                        if (value > 1)value=1;
                    }
                    _output.SetValue(k, value);

                } else {
                    _output.SetValue(k,_atlas.GetValue(k));
                }
            }

            if (denominator <= 0) {
                int x,y,z;
                _input.IndexToVoxel(i, x, y, z);
                std::cerr<<"Division by 0 while computing probabilities at voxel "<<x<<","<<y<<","<<z<<std::endl;
            }
        } else {
            for (k = 0; k < _number_of_tissues ; k++) {
                _output.SetValue(k, 0);
            }
        }
        if(_postpen)pptr++;
        ptr++;
        pm++;
        _atlas.Next();
        _output.Next();
    }
    delete[] numerator;
    delete[] G;

}








double DrawEM::Iterate(int i)
{
    if( i == 0 )
    {
        this->EStep();
    }
    else
    {
        this->EStepMRF();
    }

    if(huipvcorr){
        csflabel=0;
        gmlabel=1;
        wmlabel=2;
        this->huiPVCorrection();
    }


    this->MStep();
    std::cout << std::endl << std::endl << "After M STEP " << std::endl << std::endl;
    Print();
    this->WStep();
    this->BStep();

    std::cout << std::endl << std::endl << "After B STEP " << std::endl << std::endl;
    Print();
    return LogLikelihood();
}

void DrawEM::SetBiasField(BiasField *biasfield)
{
    _biasfield = biasfield;
}

void DrawEM::GetBiasCorrectedImage(RealImage &image)
{
    image = _input;
    //_biascorrection.Apply(_input);
}

void DrawEM::GetBiasField(RealImage &image)
{
    BytePixel *pm = _mask.GetPointerToVoxels();
    RealPixel *ptrA = _input.GetPointerToVoxels();
    RealPixel *ptrB = _uncorrected.GetPointerToVoxels();
    RealPixel *ptrOutput = image.GetPointerToVoxels();

    for( int i = 0; i < image.GetNumberOfVoxels(); ++i )
    {
        if( *pm == 1 )
        {
            *ptrOutput = *ptrB - *ptrA;
        }
        else
        {
            *ptrOutput = 0;
        }
        ptrOutput++;
        ptrA++;
        ptrB++;
        pm++;
    }
}

bool DrawEM::isPVclass(int pvclass)
{
    if( pv_classes.find(pvclass) != pv_classes.end() ) return true;
    return false;
}


void DrawEM::ConstructSegmentationHui(IntegerImage &segmentation)
{
    int i, j, m;
    RealPixel max;

    // Initialize pointers of probability maps
    _output.First();

    // Initialize segmentation to same size as input
    segmentation = IntegerImage(_input.Attributes());
    RealPixel *ptr = _input.GetPointerToVoxels();
    int *sptr = segmentation.GetPointerToVoxels();
    BytePixel *pm = _mask.GetPointerToVoxels();
    for (i = 0; i < _number_of_voxels; i++) {
        m=0;
        if (*pm == 1) {
            max  = 0;
            for (j = 0; j < _number_of_tissues; j++) {
                if (_output.GetValue(j) > max) {
                    max  = _output.GetValue(j);
                    m = tissuelabels[j];
                    if ( _has_background && (j+1) == _number_of_tissues) m=0;
                }
            }
        }
        *sptr = m;
        sptr++;
        ptr++;
        pm++;
        _output.Next();
    }
}


void DrawEM::huiPVCorrection(bool changePosterior){
    double lambda=0.5;

    if(changePosterior)lambda=0;

    std::cout<<"Hui PV correction "<<outlabel<<csflabel<<gmlabel<<wmlabel<<std::endl;
    IntegerImage segmentation, actualsegmentation;
    GreyImage scc, csfscc, outscc, mask(_input.Attributes());
    ConstructSegmentationHui(segmentation);

    ConnectivityType connectivity = (ConnectivityType) 6;
    ConnectedComponents<GreyPixel> cc(CC_LargestFirst, connectivity);



    //find small wm components
    GreyPixel *ptr = mask.GetPointerToVoxels();
    int *sptr = segmentation.GetPointerToVoxels();
    for (int j = 0; j < _number_of_voxels; j++) {
        *ptr=(*sptr==wmlabel);
        ptr++; sptr++;
    }
    cc.Input(&mask);
    cc.Output(&scc);
    cc.Run();
    int wmcomps=cc.NumberOfComponents();

    //find small csf components
    ptr = mask.GetPointerToVoxels();
    sptr = segmentation.GetPointerToVoxels();
    for (int j = 0; j < _number_of_voxels; j++) {
        *ptr=(*sptr==csflabel);
        ptr++; sptr++;
    }
    cc.Input(&mask);
    cc.Output(&csfscc);
    cc.Run();
    int csfcomps=cc.NumberOfComponents();

    //find small out components
    ptr = mask.GetPointerToVoxels();
    sptr = segmentation.GetPointerToVoxels();
    for (int j = 0; j < _number_of_voxels; j++) {
        *ptr=(*sptr==outlabel);
        ptr++; sptr++;
    }
    cc.Input (&mask);
    cc.Output(&outscc);
    cc.Run();


    int csfneighborswm[csfcomps];
    int csfneighbors[csfcomps];
    int csfneighborslven[csfcomps];
    int csfneighborsrven[csfcomps];
    for (int i=0;i<csfcomps;i++){
        csfneighborswm[i]=0;csfneighbors[i]=0;csfneighborsrven[i]=0;csfneighborslven[i]=0;
    }


    int wmvol[wmcomps];
    for(int i=0;i<wmcomps;i++) wmvol[i]=0;




    for( int x = 0; x < _input.GetX(); ++x)
    {
        for( int y = 0; y < _input.GetY(); ++y)
        {
            for( int z = 0; z < _input.GetZ(); ++z)
            {

                int comp=scc.Get(x,y,z)-1;
                wmvol[comp]++;

                comp= csfscc.Get(x,y,z)-1;


                if (comp > 0) {
                    int lx = max(x-1,0);
                    int rx = min(x+1, _input.GetX()-1 );
                    int ly = max(y-1,0);
                    int ry = min(y+1, _input.GetY()-1 );
                    int lz = max(z-1,0);
                    int rz = min(z+1, _input.GetZ()-1 );


                    if(segmentation.Get(rx ,y,z)==wmlabel)csfneighborswm[comp]++;
                    if(segmentation.Get(lx ,y,z)==wmlabel)csfneighborswm[comp]++;
                    if(segmentation.Get(x ,ry,z)==wmlabel)csfneighborswm[comp]++;
                    if(segmentation.Get(x ,ly,z)==wmlabel)csfneighborswm[comp]++;
                    if(segmentation.Get(x ,y,rz)==wmlabel)csfneighborswm[comp]++;
                    if(segmentation.Get(x ,y,lz)==wmlabel)csfneighborswm[comp]++;

                    csfneighbors[comp]+=6;
                    if(segmentation.Get(rx ,y,z)==csflabel)csfneighbors[comp]--;
                    if(segmentation.Get(lx ,y,z)==csflabel)csfneighbors[comp]--;
                    if(segmentation.Get(x ,ry,z)==csflabel)csfneighbors[comp]--;
                    if(segmentation.Get(x ,ly,z)==csflabel)csfneighbors[comp]--;
                    if(segmentation.Get(x ,y,rz)==csflabel)csfneighbors[comp]--;
                    if(segmentation.Get(x ,y,lz)==csflabel)csfneighbors[comp]--;
                }

            }
        }
    }

    _output.First();
    for( int x = 0; x < _input.GetX(); ++x)
    {
        for( int y = 0; y < _input.GetY(); ++y)
        {
            for( int z = 0; z < _input.GetZ(); ++z)
            {

                double wmval,gmval,csfval,outval;
                double owmval,ogmval,ocsfval,ooutval;

                getHuiValues(outval,csfval,gmval,wmval,x,y,z,true);
                getHuiValues(ooutval,ocsfval,ogmval,owmval,x,y,z,false);

                bool changed=false;

                int comp= csfscc.Get(x,y,z)-1;
                if (comp > 0){
                    //if is a csf cc mainly surrounded by wm -> wm
                    if(csfneighborswm[comp]>=csfneighbors[comp]/2) {
                        wmval=wmval+(1-lambda)*csfval;
                        csfval=csfval*lambda;
                        owmval=owmval+(1-lambda)*ocsfval;
                        ocsfval=ocsfval*lambda;

                        changed=true;
                    }
                }


                comp=scc.Get(x,y,z)-1;
                if(comp>0 && wmvol[comp]<0.5*wmvol[0] ){
                    csfval=csfval+(1-lambda)*(wmval+gmval);
                    wmval=wmval*lambda;
                    gmval=gmval*lambda;
                    ocsfval=ocsfval+(1-lambda)*(owmval+ogmval);
                    owmval=owmval*lambda;
                    ogmval=ogmval*lambda;
                    changed=true;
                }

                if(changed){
                    setHuiValues(outval,csfval,gmval,wmval,x,y,z,true);
                    setHuiValues(ooutval,ocsfval,ogmval,owmval,x,y,z,false);

                }
            }
        }
    }


    ConstructSegmentationHui(segmentation);
    ConstructSegmentation(actualsegmentation);


    for( int x = 0; x < _input.GetX(); ++x)
    {
        for( int y = 0; y < _input.GetY(); ++y)
        {
            for( int z = 0; z < _input.GetZ(); ++z)
            {
                if (_mask.Get(x,y,z) == 1) {

                    if(segmentation.Get(x ,y,z)!=wmlabel && segmentation.Get(x ,y,z)!=gmlabel)continue;


                    int lx = max(x-1,0);
                    int rx = min(x+1, _input.GetX()-1 );
                    int ly = max(y-1,0);
                    int ry = min(y+1, _input.GetY()-1 );
                    int lz = max(z-1,0);
                    int rz = min(z+1, _input.GetZ()-1 );

                    int neighborscsf=0;
                    if(segmentation.Get(rx ,y,z)==csflabel)neighborscsf++;
                    if(segmentation.Get(lx ,y,z)==csflabel)neighborscsf++;
                    if(segmentation.Get(x ,ry,z)==csflabel)neighborscsf++;
                    if(segmentation.Get(x ,ly,z)==csflabel)neighborscsf++;
                    if(segmentation.Get(x ,y,rz)==csflabel)neighborscsf++;
                    if(segmentation.Get(x ,y,lz)==csflabel)neighborscsf++;

                    int neighborsgm=0;
                    if(segmentation.Get(rx ,y,z)==gmlabel)neighborsgm++;
                    if(segmentation.Get(lx ,y,z)==gmlabel)neighborsgm++;
                    if(segmentation.Get(x ,ry,z)==gmlabel)neighborsgm++;
                    if(segmentation.Get(x ,ly,z)==gmlabel)neighborsgm++;
                    if(segmentation.Get(x ,y,rz)==gmlabel)neighborsgm++;
                    if(segmentation.Get(x ,y,lz)==gmlabel)neighborsgm++;

                    int neighborsout=0;
                    if(segmentation.Get(rx ,y,z)==outlabel || _mask.Get(x,y,z)==0)neighborsout++;
                    if(segmentation.Get(lx ,y,z)==outlabel || _mask.Get(x,y,z)==0)neighborsout++;
                    if(segmentation.Get(x ,ry,z)==outlabel || _mask.Get(x,y,z)==0)neighborsout++;
                    if(segmentation.Get(x ,ly,z)==outlabel || _mask.Get(x,y,z)==0)neighborsout++;
                    if(segmentation.Get(x ,y,rz)==outlabel || _mask.Get(x,y,z)==0)neighborsout++;
                    if(segmentation.Get(x ,y,lz)==outlabel || _mask.Get(x,y,z)==0)neighborsout++;

                    double wmval,gmval,csfval,outval;
                    double owmval,ogmval,ocsfval,ooutval;


                    getHuiValues(outval,csfval,gmval,wmval,x,y,z,true);
                    getHuiValues(ooutval,ocsfval,ogmval,owmval,x,y,z,false);

                    bool changed=false;



                    if(segmentation.Get(x ,y,z)==wmlabel){
                        if(wmval==0)wmval=0.1;
                        //if is a wm voxel that touches outlier and csf -> csf
                        if( neighborsout>0 && neighborscsf>0){
                            csfval=csfval+(1-lambda)*(wmval+gmval);
                            wmval=wmval*lambda;
                            gmval=gmval*lambda;

                            ocsfval=ocsfval+(1-lambda)*(owmval+ogmval);
                            owmval=owmval*lambda;
                            ogmval=ogmval*lambda;

                            changed=true;
                        }else if((neighborscsf>0 || neighborsout>0) && neighborsgm>0){
                            //if is a wm voxel that touches (outlier or csf) and gm -> csf, gm

                            double sum=gmval+csfval;
                            if(sum!=0){
                                gmval=gmval+(1-lambda)*wmval * (gmval/sum);
                                csfval=csfval+(1-lambda)*wmval * (csfval/sum);
                            }else{
                                gmval=gmval+(1-lambda)*0.5 * wmval;
                                csfval=csfval+(1-lambda)*0.5 * wmval;
                            }
                            wmval=wmval*lambda;


                            double osum=ogmval+ocsfval;
                            if(osum!=0){
                                ogmval=ogmval+(1-lambda)*owmval * (ogmval/osum);
                                ocsfval=ocsfval+(1-lambda)*owmval * (ocsfval/osum);
                            }else{
                                ogmval=ogmval+(1-lambda)*0.5 * owmval;
                                ocsfval=ocsfval+(1-lambda)*0.5 * owmval;
                            }
                            owmval=owmval*lambda;

                            changed=true;
                        }

                    }else{
                        if(gmval==0)gmval=0.1;
                        if(neighborsout>0){
                            if(neighborscsf>0){
                                //if is a gm voxel that touches outlier and csf -> csf
                                csfval=csfval+(1-lambda)*(wmval+gmval);
                                ocsfval=ocsfval+(1-lambda)*(owmval+ogmval);
                                changed=true;
                            }
                            else if(neighborsgm==0){
                                //if is a gm voxel that touches outlier and does not have gm neighbors -> out
                                outval=outval+(1-lambda)*(wmval+gmval);
                                ooutval=ooutval+(1-lambda)*(owmval+ogmval);
                                changed=true;
                            }
                            if(changed){
                                wmval=wmval*lambda;
                                gmval=gmval*lambda;
                                owmval=owmval*lambda;
                                ogmval=ogmval*lambda;

                            }
                        }
                    }


                    if(changed){
                        setHuiValues(outval,csfval,gmval,wmval,x,y,z,true);
                        setHuiValues(ooutval,ocsfval,ogmval,owmval,x,y,z,false);
                    }

                }
            }
        }
    }


}





void DrawEM::getHuiValues(double &outval,double &csfval,double &gmval,double &wmval,int x,int y,int z,bool atlas){
    double vals[5];
    for(int i=0;i<5;i++)vals[i]=0;

    for(int j=0;j<_number_of_tissues;j++){
        if(tissuelabels[j]==0 )continue;
        if(atlas)vals[tissuelabels[j]]+=_atlas.GetValue(x,y,z,j);
        else     vals[tissuelabels[j]]+=_output.GetValue(x,y,z,j);
    }

    outval=vals[outlabel];
    csfval=vals[csflabel];
    gmval=vals[gmlabel];
    wmval=vals[wmlabel];
}

void DrawEM::setHuiValues(double &outval,double &csfval,double &gmval,double &wmval,int x,int y,int z,bool atlas){
    double vals[5],newvals[5],hm[5];
    for(int i=0;i<5;i++){
        vals[i]=0;hm[i]=0;
    }

    newvals[outlabel]=outval;
    newvals[csflabel]=csfval;
    newvals[gmlabel]=gmval;
    newvals[wmlabel]=wmval;

    for(int j=0;j<_number_of_tissues;j++){
        if(tissuelabels[j]==0 )continue;
        if(atlas)vals[tissuelabels[j]]+=_atlas.GetValue(x,y,z,j);
        else     vals[tissuelabels[j]]+=_output.GetValue(x,y,z,j);
        hm[tissuelabels[j]]++;
    }

    double val,newval,part;
    for(int j=0;j<_number_of_tissues;j++){
        if(tissuelabels[j]==0 )continue;

        if(atlas)val=_atlas.GetValue(x,y,z,j);
        else val=_output.GetValue(x,y,z,j);

        if(vals[tissuelabels[j]]==0)
            part=1/hm[tissuelabels[j]];
        else
            part=val/vals[tissuelabels[j]];
        newval=part*newvals[tissuelabels[j]];


        if(atlas)
            _atlas.SetValue(x,y,z,j,newval);
        else
            _output.SetValue(x,y,z,j, newval);
    }
}




double DrawEM::getMRFInterEnergy(int index, int tissue)
{
    if( _connectivity.Rows() == 1 || !intermrf)
    {
        return 1.0;
    }

    double expo = -1.0* mrfweight;

    double energy = .0;
    int x,y,z;
    _input.IndexToVoxel(index, x, y, z);

    expo *= betainter;


    double weight = 0;

    for( int k = 0; k < _number_of_tissues; k++)
    {
        double temp = 0;
        weight = 0;
        double conn = _connectivity.Get(k, tissue);

        weight=conn;
        if (conn==0) continue;

        temp += (_MRF_inter[k])->Get( x ,y,z);
        energy += temp * weight;
    }


    expo *= energy;

    return exp(expo);
}


template DrawEM::DrawEM(int, RealImage **, RealImage *);
template DrawEM::DrawEM(int, HashRealImage **, HashRealImage *);
template DrawEM::DrawEM(int, RealImage **);
template DrawEM::DrawEM(int, HashRealImage **);
template DrawEM::DrawEM(int, RealImage **, RealImage **);
template DrawEM::DrawEM(int, HashRealImage **, HashRealImage **);

}
