function d = density_cloud(pts)
%DENSITY
% d = density(coords) 
% Returns density equal to the sum of distances to the two nearest neighbors
% in a fixed point cloud, divided by 20.
% Takes matrices of size dim x N, where N is the number of vectors. 
m =    [   -0.7846   -0.3702   -0.5450;
           -0.5553   -0.2670   -0.7755;
           -0.4057   -0.7614   -0.0445;
           -0.6212   -0.0156   -0.7291;
            0.0207   -0.1215    0.3273;
            0.6826   -0.1796   -0.2404;
            0.0775    0.4723   -0.6585;
            0.2982   -0.2584    0.2813;
            0.5155   -0.3749   -0.3344;
           -0.5983   -0.3355    0.2421;
           -0.0117   -0.7436   -0.9377;
            0.9439   -0.7930   -0.6228;
            0.5615    0.3719    0.9117;
           -0.0426   -0.1283   -0.3705;
            0.5054   -0.9873    0.0682;
            0.8931    0.4838   -0.1078;
           -0.6272    0.8160    0.0185;
           -0.9291    0.4349    0.3433;
            0.9273    0.0046    0.4035;
           -0.1073   -0.7305   -0.7849;
           -0.1722   -0.4491    0.3682;
            0.0766    0.1684   -0.0575;
           -0.8580   -0.0663   -0.6879;
            0.3306    0.0290   -0.8627;
           -0.3240   -0.5483    0.8015;
            0.3277   -0.0300   -0.2903;
            0.0691    0.1532   -0.9564;
           -0.8199   -0.6919   -0.5147;
            0.2613   -0.0957    0.5349;
           -0.0232   -0.8705    0.4106;
            0.3385   -0.5726    0.1589;
           -0.7733    0.9008   -0.1762;
           -0.4221    0.9833    0.5486;
           -0.1023   -0.5968   -0.2547;
           -0.6769    0.7268   -0.2219;
           -0.9785    0.2244   -0.9543;
           -0.2894   -0.0456   -0.8039;
            0.8298    0.6202    0.2574;
            0.5802   -0.2560    0.8264;
            0.9817    0.7439   -0.3070;
            0.1547    0.8400    0.6856;
           -0.1122   -0.0766    0.2103;
            0.3574   -0.0359    0.7831;
            0.4431   -0.6322    0.3271;
            0.1273    0.1940    0.7562;
            0.4704    0.7640   -0.3234;
           -0.5414   -0.1907   -0.3315;
            0.9653   -0.6201    0.0209;
           -0.2277   -0.0018    0.2840;
            0.2916    0.1014    0.4479;
            0.7408   -0.4255   -0.6174;
           -0.2693    0.2457   -0.8210;
           -0.9285    0.2184   -0.0930;
           -0.7380    0.2969   -0.9009;
           -0.1439    0.8777    0.0845;
           -0.4863   -0.1949   -0.6704;
           -0.8977    0.6233    0.5685;
           -0.7452    0.7760   -0.4784;
            0.3869   -0.4791   -0.5420;
            0.7611   -0.9320    0.6767;
           -0.2800   -0.5526    0.1172;
           -0.6874   -0.1671    0.1516;
           -0.5130   -0.2423   -0.3423;
            0.9650   -0.1292   -0.0188;
           -0.5203    0.6470    0.1979;
           -0.2514    0.9558    0.4971;
           -0.3916    0.9818   -0.7191;
           -0.0162    0.6169   -0.7067;
           -0.5468    0.5129    0.9624;
           -0.3373   -0.5770   -0.3465;
            0.3323   -0.6967   -0.9373;
           -0.7580   -0.1619    0.5770;
            0.7076   -0.5769    0.0347;
            0.9943   -0.2420    0.2465;
            0.2035   -0.7617    0.0672;
           -0.9802   -0.9929   -0.2356;
           -0.9510   -0.2822    0.8173;
           -0.9616    0.5846   -0.4899;
            0.5911    0.5657    0.5516;
           -0.0213    0.6937    0.4937;
            0.8865   -0.3867    0.7670;
            0.8857   -0.9067    0.2069;
            0.9213    0.6117    0.6429;
           -0.8500   -0.6126    0.3148;
           -0.2203   -0.8267    0.4757;
           -0.8808   -0.4385   -0.1643;
            0.5450   -0.9680   -0.0929;
            0.6201    0.4168   -0.6268;
            0.9649    0.5441   -0.1712;
           -0.4107   -0.8789    0.2148;
            0.3671   -0.6443    0.5199;
           -0.2822    0.9463    0.1034;
           -0.3885   -0.9474    0.7454;
            0.9141   -0.8353    0.3859;
            0.7759   -0.4758   -0.0683;
           -0.5237   -0.5902   -0.2464;
           -0.0051    0.6405    0.4367;
            0.9183    0.7985   -0.4037;
           -0.0250   -0.3798   -0.4489;
            0.7465    0.4446    0.0038] ;
[~, dis] = knnsearch(m, pts','k',2);
d = sum(dis,2)' / 20;