/**
 *
 * @author Timothy B. Hayward
 */

package extended_kinematic_fitters;

import org.jlab.clas.physics.GenericKinematicFitter;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;

import org.jlab.clas.physics.*;

//import org.apache.commons.math3.*;

public class analysis_fitter extends GenericKinematicFitter {

    protected final Double mybeam;

    public analysis_fitter(double beam) {
        super(beam);
        mybeam = beam;
    } 
    
    ////////////////////////////////////////////////////////////////////////////////////////////////  
    // General cuts to test the validity of the event and particle
    public boolean banks_test(DataEvent event) {
        return event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter") && event.hasBank("REC::Track") && 
            event.hasBank("REC::Traj") && event.hasBank("REC::Cherenkov") &&
            event.hasBank("RUN::config");
    }
    
    public boolean forward_detector_cut(int current_Part, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", current_Part);
//        System.out.println(status);
        return (Math.abs(status)<4000 || Math.abs(status)>4999) && Math.abs(status)>1999;
    }
    
    public boolean highest_e_in_fd_cut(int p_max_index, HipoDataBank rec_Bank) {
        int status = rec_Bank.getInt("status", p_max_index);
        return status<-1999;
    }
    
    public boolean theta_cut(int current_Part, HipoDataBank rec_Bank) {
        float px = rec_Bank.getFloat("px", current_Part);
        float py = rec_Bank.getFloat("py", current_Part);
        float pz = rec_Bank.getFloat("pz", current_Part);
        double r = Math.pow(px*px + py*py + pz*pz, 0.5);
        double theta = (180/Math.PI)*Math.acos(pz/r);
        
        return 6<theta && theta<30;
//        return 15<theta && theta<25;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Fiducial volume cuts on the pcal and drift chambers
    
    public boolean pcal_fiducial_cut(int current_Part, HipoDataBank cal_Bank) {
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            if (cal_Bank.getInt("pindex", current_Row)==current_Part && 
                    cal_Bank.getInt("layer", current_Row)==1) {
                float lv = cal_Bank.getFloat("lv", current_Row);
                float lw = cal_Bank.getFloat("lw", current_Row);
                return lv > 9 && lw > 9; // "loose cut, change to 14 and 19 for medium and tight"
                
            }
        }
        return false;
    }

    public boolean dc_fiducial_cut(int current_Part, HipoDataBank rec_Bank, 
    HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        // trajectory crosses  (layers: 6 = DC region 1 start,  18 = DC region 2 center,  
        // 36 = DC region 3 end)
        
        
        int event_number = run_Bank.getInt("event", 0);
        //fitted values
    double[][][][] maxparams_had = {
        {{{-51.2054, 38.3199, -2.34912, 0.0195708},{-35.124, 26.5405, -1.09976, 0.0054436},{-38.4667, 28.8007, -1.29647, 0.00723946}},
{{-37.5758, 26.6792, -1.00373, 0.00425542},{-45.5585, 35.2513, -2.21029, 0.0196454},{-40.6878, 30.355, -1.42152, 0.00817034}},
{{-38.3878, 29.6034, -1.68592, 0.0150766},{-36.8921, 27.4735, -1.14582, 0.00554924},{-23.5149, 18.0147, -0.38597, 2.87092e-09}},
{{-48.2748, 35.4465, -1.92023, 0.0136412},{-32.9275, 24.1143, -0.715458, 5.40587e-05},{-62.2819, 53.7417, -5.16183, 0.0613377}},
{{-39.9916, 30.4254, -1.72002, 0.015037},{-20.4886, 16.3741, -0.310525, 2.68397e-08},{-38.2334, 31.4077, -2.12755, 0.021869}},
{{-62.9673, 48.747, -3.71315, 0.036982},{-50.8588, 37.2435, -1.98066, 0.0126239},{-72.3186, 58.7791, -5.21921, 0.0573224}}},
{{{-7.47342e-07, 12.8291, -0.819744, 0.00818184},{-3.46906, 12.9438, -0.594252, 0.00438008},{-6.57983, 11.8439, -0.240227, 3.27439e-10}},
{{-2.23025e-08, 12.8274, -0.819395, 0.00817682},{-4.95047, 13.5494, -0.600411, 0.00406046},{-8.411, 12.7605, -0.278882, 2.14459e-10}},
{{-3.62667e-08, 13.2965, -0.904711, 0.00919042},{-9.89194, 16.9776, -0.934777, 0.00761783},{-5.88727, 11.4158, -0.218181, 5.04197e-11}},
{{-1.51214e-08, 12.7606, -0.797213, 0.00785312},{-5.63162, 14.3101, -0.702664, 0.00528812},{-4.77062, 10.9184, -0.200027, 2.67841e-13}},
{{-1.81438e-08, 12.9692, -0.848261, 0.00859668},{-8.30238, 16.2372, -0.899398, 0.00745394},{-7.58827, 12.4473, -0.270602, 2.53855e-10}},
{{-4.29779e-08, 13.121, -0.861715, 0.0085118},{-7.36773, 15.334, -0.776242, 0.00585805},{-6.54892, 12.0405, -0.259209, 3.73643e-06}}},
{{{-2.68279e-07, 12.99, -0.846226, 0.00845788},{-14.6317, 19.3874, -1.09244, 0.00899541},{-38.1915, 29.8688, -1.59229, 0.0120089}},
{{-0.996514, 13.9379, -0.964686, 0.00982941},{-15.9613, 20.2461, -1.16106, 0.00955431},{-35.9455, 29.0996, -1.586, 0.0122175}},
{{-1.14284e-07, 13.6015, -0.966952, 0.0101523},{-15.5288, 20.3045, -1.20523, 0.0102808},{-34.2682, 26.4216, -1.20609, 0.0078434}},
{{-1.70075e-08, 13.0005, -0.832325, 0.00817159},{-7.66776, 15.4526, -0.779727, 0.00585967},{-26.8035, 23.9995, -1.2322, 0.00942061}},
{{-9.53804e-10, 13.2563, -0.898206, 0.00917629},{-6.85083, 14.8485, -0.722803, 0.0053221},{-39.3606, 31.5412, -1.83015, 0.0148302}},
{{-7.66835e-07, 13.937, -1.05153, 0.0118223},{-9.7913, 16.925, -0.913158, 0.00712552},{-27.722, 23.9412, -1.1314, 0.00761088}}},
{{{-22.1832, 20.4134, -0.764848, 0.00310923},{-31.0844, 28.2369, -1.715, 0.0145145},{-9.52175, 18.7932, -1.38896, 0.0150233}},
{{-21.5849, 20.2457, -0.762109, 0.00305359},{-19.5601, 21.5945, -1.18955, 0.00939109},{-1.57084, 13.3989, -0.823161, 0.00795227}},
{{-16.052, 16.6264, -0.444308, 2.82701e-06},{-13.8291, 18.6541, -1.01549, 0.00825776},{-1.92223e-05, 13.0305, -0.881089, 0.00925281}},
{{-19.821, 18.4301, -0.516168, 2.17199e-10},{-30.6295, 28.0989, -1.71897, 0.0146585},{-9.23709, 17.1589, -1.03955, 0.00943673}},
{{-16.1795, 16.7121, -0.448883, 1.53774e-11},{-23.6418, 24.5748, -1.48652, 0.01254},{-4.2626e-09, 12.899, -0.845374, 0.00872171}},
{{-9.74791, 15.0287, -0.531727, 0.00192371},{-41.0848, 33.1802, -1.97671, 0.0158148},{-4.12428, 14.3361, -0.820483, 0.00725632}}},
{{{-106.938, 86.1515, -8.33159, 0.101803},{-147.879, 123.116, -13.5824, 0.189579},{-40.461, 32.3475, -1.89454, 0.0152611}},
{{-116.247, 96.2163, -9.93814, 0.130211},{-148.54, 117.766, -12.0624, 0.156313},{-21.2546, 20.4237, -0.895905, 0.00581147}},
{{-1.80778e-07, 13.0131, -0.835058, 0.00806212},{-142.768, 110.58, -10.8081, 0.133187},{-68.8917, 52.0911, -3.89921, 0.0389652}},
{{-120.479, 92.6075, -8.56746, 0.0999391},{-147.549, 115.807, -11.6707, 0.149171},{-68.1786, 52.4364, -4.06328, 0.042723}},
{{-123.541, 97.8256, -9.57399, 0.118002},{-149.955, 123.896, -13.5475, 0.187947},{-17.0127, 16.3323, -0.393286, 1.67366e-13}},
{{-98.3311, 81.6497, -8.10227, 0.102672},{-149.833, 127.44, -14.5373, 0.209948},{-0.74045, 9.36883, -0.152603, 6.37437e-05}}},
{{{-142.357, 114.799, -11.8528, 0.152895},{-149.875, 111.604, -10.4642, 0.125095},{-149.501, 121.578, -13.0379, 0.177271}},
{{-29.9249, 22.0873, -0.60085, 1.3473e-10},{-42.0584, 33.3444, -1.99565, 0.016557},{-150, 127.502, -14.6523, 0.21559}},
{{-118.755, 95.587, -9.52135, 0.120561},{-144.344, 109.666, -10.4534, 0.126138},{-149.437, 93.2368, -6.50292, 0.0579289}},
{{-149.814, 110.828, -10.2785, 0.122489},{-148.167, 114.308, -11.2695, 0.140999},{-149.539, 112.153, -10.6506, 0.128573}},
{{-33.521, 26.8939, -1.29466, 0.00851674},{-141.738, 107.08, -10.1005, 0.120796},{-150, 108.507, -9.8741, 0.116558}},
{{-26.8311, 20.6856, -0.549123, 5.34572e-11},{-148.673, 110.129, -10.1547, 0.118973},{-120.868, 101.527, -10.6726, 0.13942}}}};

    double[][][][] minparams_had = {
        {{{38.6069, -28.188, 1.18288, -0.00589988},{26.3209, -21.6986, 0.768368, -0.00281998},{15.8707, -11.8368, 3.22688e-09, -8.36654e-10}},
{{44.6032, -32.8751, 1.65153, -0.0105575},{51.7598, -35.5217, 1.29448, -2.11513e-06},{32.9204, -26.2522, 1.46865, -0.0139537},},
{{32.1412, -24.2383, 1.00154, -0.0064333},{31.3784, -23.4733, 0.774187, -0.00203468},{53.8936, -44.9217, 3.8766, -0.0442173}},
{{43.416, -32.3319, 1.78081, -0.0143004},{50.0314, -38.5805, 2.56055, -0.0235214},{25.8309, -21.3241, 1.01076, -0.00944255}},
{{50.4741, -40.1723, 3.0563, -0.0325218},{61.6178, -50.3165, 4.34833, -0.0483999},{12.7766, -11.4712, 0.0511044, -1.74577e-10}},
{{43.2153, -31.5115, 1.4834, -0.00863028},{69.3381, -51.2282, 3.39516, -0.0272937},{120.36, -85.4755, 7.45064, -0.0836405}}},
{{{0.923008, -13.7147, 0.895123, -0.00853805},{7.6279, -15.4483, 0.798758, -0.00633026},{5.62644, -11.4345, 0.226312, -1.25172e-09}},
{{2.50458, -14.7412, 1.01929, -0.0102917},{9.0023, -16.3615, 0.894467, -0.00751183},{3.69192, -10.4985, 0.188359, -3.1185e-10}},
{{1.5978, -13.7835, 0.890884, -0.00868481},{13.4977, -19.2788, 1.15784, -0.0101226},{3.82792, -10.4417, 0.179639, -2.24507e-10}},
{{1.97078e-05, -13.0168, 0.817426, -0.00762709},{6.07517, -14.4662, 0.696387, -0.00504973},{6.98517, -12.2129, 0.264062, -2.7665e-11}},
{{2.39018e-09, -13.0476, 0.838641, -0.00809971},{4.15408, -14.0986, 0.777823, -0.00678748},{3.77333, -10.6691, 0.201849, -4.83515e-10}},
{{0.00103141, -13.1638, 0.857402, -0.00828706},{7.16766, -15.2964, 0.811382, -0.00673444},{3.39768, -10.2477, 0.172417, -3.86335e-10}}},
{{{1.59369e-06, -13.8294, 0.990918, -0.0103128},{20.1273, -23.853, 1.58449, -0.0145959},{40.8152, -32.8944, 2.00731, -0.0171007}},
{{1.4334, -14.5452, 1.04379, -0.0106791},{19.9242, -23.3894, 1.5036, -0.0134429},{45.1348, -34.9897, 2.11238, -0.0175613}},
{{4.48276e-06, -12.6688, 0.757818, -0.006981},{10.2525, -16.9056, 0.909637, -0.00739798},{33.2958, -27.7763, 1.53467, -0.0123488}},
{{3.817e-06, -13.2285, 0.856439, -0.0081744},{12.5356, -19.0801, 1.1686, -0.0102758},{37.3388, -29.7344, 1.64296, -0.0130658}},
{{3.64842e-07, -14.1631, 1.0771, -0.0118569},{9.85442, -17.8198, 1.12641, -0.010627},{34.7, -28.5335, 1.57226, -0.0124004}},
{{0.828721, -13.6429, 0.895665, -0.00866683},{10.8176, -18.0919, 1.11147, -0.010183},{29.9288, -24.3389, 1.08973, -0.00703934}}},
{{{15.8302, -16.9632, 0.53561, -0.00136216},{32.8002, -29.2569, 1.79783, -0.015324},{1.98393, -13.0099, 0.70788, -0.00615153}},
{{16.0367, -16.5901, 0.470678, -0.000728065},{32.4005, -29.7403, 1.92286, -0.0171968},{2.39707, -13.6612, 0.816883, -0.00770837}},
{{22.0623, -21.6319, 1.02811, -0.00680893},{32.7467, -29.6099, 1.87839, -0.0164223},{1.19902e-08, -12.972, 0.863127, -0.00884759}},
{{21.5883, -21.198, 0.957819, -0.00575361},{25.7387, -25.4963, 1.5428, -0.0131855},{6.06479, -16.6311, 1.16092, -0.0117194}},
{{19.6915, -19.1751, 0.704086, -0.00288768},{28.6596, -27.3351, 1.70309, -0.0148193},{5.30096e-08, -11.8562, 0.621373, -0.00541869}},
{{20.6594, -19.8704, 0.786033, -0.00394155},{20.7612, -22.3774, 1.27116, -0.0104109},{2.56196, -14.4159, 0.98009, -0.0100214}}},
{{{14.2631, -23.4009, 1.97498, -0.0231565},{4.79173, -12.3009, 0.326886, -1.01998e-11},{116.734, -88.2089, 8.05103, -0.0957363}},
{{5.23149, -16.0375, 1.0661, -0.0101511},{110.165, -85.9607, 7.97299, -0.0942078},{43.5095, -33.7305, 1.93256, -0.0149425}},
{{1.85579e-09, -11.9228, 0.619166, -0.00523146},{7.01255, -15.5132, 0.84916, -0.00706963},{149.945, -108.603, 9.82162, -0.114667}},
{{140.48, -115.163, 12.1689, -0.16217},{133.49, -104.362, 10.1931, -0.125813},{114.099, -83.4723, 7.13202, -0.0793527}},
{{18.6833, -24.9619, 1.82946, -0.0167602},{11.0804, -17.1484, 0.877833, -0.00652953},{15.5084, -15.7594, 0.378209, -2.63283e-06}},
{{4.64372e-08, -8.1395, 0.0139504, -9.40293e-06},{74.1514, -64.5468, 6.22895, -0.0765922},{29.0616, -24.2817, 1.13928, -0.00789879}}},
{{{150, -119.236, 12.4315, -0.16672},{17.764, -15.9517, 0.338355, -9.39693e-11},{149.387, -101.062, 8.19036, -0.0851116}},
{{130.018, -97.9952, 9.09323, -0.109458},{147.831, -109.347, 10.1054, -0.119381},{149.958, -109.621, 10.0894, -0.119333}},
{{149.744, -119.521, 12.5308, -0.168366},{105.835, -82.9887, 7.74468, -0.0929233},{135.964, -117.822, 13.5665, -0.197518}},
{{150, -119.218, 12.431, -0.167282},{43.2244, -35.0239, 2.19659, -0.0185046},{150, -115.904, 11.6807, -0.151716}},
{{149.062, -108.704, 9.81712, -0.113862},{144.909, -115.147, 11.8228, -0.153516},{150, -115.103, 11.5315, -0.149515}},
{{149.245, -111.727, 10.5343, -0.126258},{127.9, -98.607, 9.4799, -0.117458},{135.306, -106.942, 10.7421, -0.136543}}}};

   
   double[][][][] maxparams_elec =
        {{{{-14.563, 0.60032},{-19.6768, 0.58729},{-22.2531, 0.544896}},
{{-12.7486, 0.587631},{-18.8093, 0.571584},{-19.077, 0.519895}},
{{-11.3481, 0.536385},{-18.8912, 0.58099},{-18.8584, 0.515956}},
{{-10.7248, 0.52678},{-18.2058, 0.559429},{-22.0058, 0.53808}},
{{-16.9644, 0.688637},{-17.1012, 0.543961},{-21.3974, 0.495489}},
{{-13.4454, 0.594051},{-19.4173, 0.58875},{-22.8771, 0.558029}}},
{{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
{{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
{{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
{{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
{{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}}};


    double[][][][] minparams_elec =
        {{{{12.2692, -0.583057},{17.6233, -0.605722},{19.7018, -0.518429}},
{{12.1191, -0.582662},{16.8692, -0.56719},{20.9153, -0.534871}},
{{11.4562, -0.53549},{19.3201, -0.590815},{20.1025, -0.511234}},
{{13.202, -0.563346},{20.3542, -0.575843},{23.6495, -0.54525}},
{{12.0907, -0.547413},{17.1319, -0.537551},{17.861, -0.493782}},
{{13.2856, -0.594915},{18.5707, -0.597428},{21.6804, -0.552287}}},
{{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
{{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
{{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
{{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}},
{{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}}};
;


        boolean track_success = true; 
        for (int current_Row = 0; current_Row < traj_Bank.rows(); current_Row++) {
            if (!track_success) { continue; }
            // loop over all entries in the trajectory bank
            if (current_Part == traj_Bank.getInt("pindex", current_Row)) {
                // require that the particle examined corresponds to this track
                int region = -1;
                if (traj_Bank.getInt("layer", current_Row) == 6) {
                    region = 0;
                } else if (traj_Bank.getInt("layer", current_Row) == 18) {
                    region = 1;
                } else if (traj_Bank.getInt("layer", current_Row) == 36) {
                    region = 2;
                } // note these are CLAS region - 1 because of java indexing rules

                double theta_DCr = 5000;
                double phi_DCr_raw = 5000;
                if (region == 0 || region == 1 || region == 2) {
                    float cx = traj_Bank.getFloat("cx", current_Row);
                    float cy = traj_Bank.getFloat("cy", current_Row);
                    float x = traj_Bank.getFloat("x", current_Row);
                    float y = traj_Bank.getFloat("y", current_Row);
                    float z = traj_Bank.getFloat("z", current_Row);
                    float cz = traj_Bank.getFloat("cz", current_Row);
                    double r = Math.sqrt(x*x + y*y + z*z);
                    theta_DCr = (180/Math.PI) * Math.acos(z/r);
                    phi_DCr_raw = (180/Math.PI) * Math.atan2(y/r,x/r);
                    
                    int sector = 0;
                    if(phi_DCr_raw < 30 && phi_DCr_raw >= -30){        sector = 0;} // note this is CLAS sector - 1
                    else if(phi_DCr_raw < 90 && phi_DCr_raw >= 30){    sector = 1;} // because of Java indexing reasons
                    else if(phi_DCr_raw < 150 && phi_DCr_raw >= 90){   sector = 2;}
                    else if(phi_DCr_raw >= 150 || phi_DCr_raw < -150){ sector = 3;}
                    else if(phi_DCr_raw < -90 && phi_DCr_raw >= -150){ sector = 4;}
                    else if(phi_DCr_raw < -30 && phi_DCr_raw >= -90){  sector = 5;}
                    double phi_DCr = 5000;
                    if (sector == 0) {phi_DCr = phi_DCr_raw;}
                    if (sector == 1) {phi_DCr = phi_DCr_raw - 60;}
                    if (sector == 2) {phi_DCr = phi_DCr_raw - 120;}
                    if (sector == 3 && phi_DCr_raw > 0) {phi_DCr = phi_DCr_raw - 180;}
                    if (sector == 3 && phi_DCr_raw < 0) {phi_DCr = phi_DCr_raw + 180;}
                    if (sector == 4) {phi_DCr = phi_DCr_raw + 120;}
                    if (sector == 5) {phi_DCr = phi_DCr_raw + 60;}
                    
                    
//                    System.out.println(phi_DCr_raw+" "+phi_DCr);

                    int pid = rec_Bank.getInt("pid", current_Part);
                    switch (pid) { // these are to set the index in the paramter lists
                    case 11:
                        pid = 0;
                        break;

                    case 2212:
                        pid = 1;
                        break;

                    case 211:
                        pid = 2;
                        break;

                    case -211:
                        pid = 3;
                        break;

                    case 321:
                        pid = 4;
                        break;

                    case -321:
                        pid = 5;
                        break;
                    }
                    
//                    if (pid==0 && sector!=0) {
//                        return false;
//                    }

                    int pdg_pid = rec_Bank.getInt("pid", current_Part);
                    if (rec_Bank.getInt("pid", current_Part) != 11) { // hadrons
                        double calc_phi_min = minparams_had[pid][sector][region][0] + 
                            minparams_had[pid][sector][region][1]*Math.log(theta_DCr) + 
                            minparams_had[pid][sector][region][2]*theta_DCr +
                            minparams_had[pid][sector][region][3]*theta_DCr*theta_DCr;

                        double calc_phi_max = maxparams_had[pid][sector][region][0] + 
                            maxparams_had[pid][sector][region][1]*Math.log(theta_DCr) + 
                            maxparams_had[pid][sector][region][2]*theta_DCr +
                            maxparams_had[pid][sector][region][3]*theta_DCr*theta_DCr; 

                        
                        track_success = phi_DCr > calc_phi_min && phi_DCr < calc_phi_max;
//                        if (track_success == false) {
//                        System.out.println("event: "+event_number+", pid: "+pdg_pid+", sector: "+(sector+1)+", region: "+(region+1)+", theta_DCr: "+theta_DCr+", calc_phi_min: "+calc_phi_min+", phi_DCr: "+phi_DCr+", calc_phi_max: "+calc_phi_max+", track success: "+track_success);
//                        }
                    } else { // electrons
                        double x_New = 10000;
                        double y_New = 10000;
                        switch (sector) {
                            case 0:
                                x_New = x;
                                y_New = y;
                                break;
                            case 1:
                                x_New = x*Math.cos(-60*Math.PI/180) - y*Math.sin(-60*Math.PI/180);
                                y_New = x*Math.sin(-60*Math.PI/180) + y*Math.cos(-60*Math.PI/180);      
                                break;
                            case 2:
                                x_New = x*Math.cos(-120*Math.PI/180) - y*Math.sin(-120*Math.PI/180);
                                y_New = x*Math.sin(-120*Math.PI/180) + y*Math.cos(-120*Math.PI/180);      
                                break;
                            case 3:
                                x_New = x*Math.cos(-180*Math.PI/180) - y*Math.sin(-180*Math.PI/180);
                                y_New = x*Math.sin(-180*Math.PI/180) + y*Math.cos(-180*Math.PI/180);      
                                break;
                            case 4:
                                x_New = x*Math.cos(120*Math.PI/180) - y*Math.sin(120*Math.PI/180);
                                y_New = x*Math.sin(120*Math.PI/180) + y*Math.cos(120*Math.PI/180);      
                                break;
                            case 5:
                                x_New = x*Math.cos(60*Math.PI/180) - y*Math.sin(60*Math.PI/180);
                                y_New = x*Math.sin(60*Math.PI/180) + y*Math.cos(60*Math.PI/180);      
                                break;
                        }
                        double calc_min = minparams_elec[pid][sector][region][0]+minparams_elec[pid][sector][region][1]*x_New;
                        double calc_max = maxparams_elec[pid][sector][region][0]+maxparams_elec[pid][sector][region][1]*x_New;
                        
                        track_success = y_New > calc_min && y_New < calc_max;
//                        System.out.println("event: "+event_number+", pid: "+pdg_pid+", sector: "+(sector+1)+", region: "+(region+1)+", x: "+x+", y: "+y+", calc_min: "+calc_min+", y_New: "+y_New+", calc_max: "+calc_max+", track success: "+track_success);
                    }
                    
                
                }
            }

        }
        return track_success;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // PID enhancements
    
    public boolean nphe_cut(int current_Part, HipoDataBank cc_Bank) {
        for (int current_Row = 0; current_Row < cc_Bank.rows(); current_Row++) {
            if (cc_Bank.getInt("pindex", current_Row)==current_Part) {
                return cc_Bank.getFloat("nphe", current_Row) > 2;
            }
        }
        return false; 
    }
    
    public boolean calorimeter_energy_cut(int current_Part, HipoDataBank cal_Bank) {
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            if (cal_Bank.getInt("pindex", current_Row)==current_Part && 
                    cal_Bank.getInt("layer", current_Row)==1) {
                return cal_Bank.getFloat("energy", current_Row) > 0.07;
            }
        }
        return false;
    }
    
    public boolean calorimeter_sampling_fraction_cut(int current_Part, double p, HipoDataBank cal_Bank) {
        double scale = 3.5; // how many std away from mean to cut on
        int sector = -1;
        double cal_energy = 0;
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            if (cal_Bank.getInt("pindex", current_Row)==current_Part)  {
                sector = cal_Bank.getInt("sector", current_Row) - 1; // subtract one to start at index of 0
                cal_energy+= cal_Bank.getFloat("energy", current_Row);
            }
        }
        
        double e_cal_sampl_mu[][] = {{0.2531, 0.2550, 0.2514, 0.2494, 0.2528, 0.2521},
            {-0.6502, -0.7472, -0.7674, -0.4913, -0.3988, -0.703},
            {4.939, 5.350, 5.102, 6.440, 6.149, 4.957}};
        
        double e_cal_sampl_sigma[][] = {{0.002726, 0.004157, 0.00522, 0.005398, 0.008453, 0.006553}, 
            {1.062, 0.859, 0.5564, 0.6576, 0.3242, 0.4423}, 
            {-4.089, -3.318, -2.078, -2.565, -0.8223, -1.274}};
        
        double mean = e_cal_sampl_mu[0][sector]+(e_cal_sampl_mu[1][sector]/1000)*(p-e_cal_sampl_mu[2][sector])*
                (p-e_cal_sampl_mu[2][sector]);
        
        double std = e_cal_sampl_sigma[0][sector] + e_cal_sampl_sigma[1][sector] / 
                (10 * (p-e_cal_sampl_sigma[2][sector]));
        
        // sampling fraction is cal_energy/p
        return ((cal_energy/p) > (mean-scale*std)) && ((cal_energy/p) < (mean+scale*std));
    }
    
    public boolean calorimeter_diagonal_cut(int current_Part, double p, HipoDataBank cal_Bank) {
        if (p < 4.5) {
            return true; // only apply diagonal cut above 4.5 GeV
        }
        double pcal_plus_ecal_inner = 0;
        for (int current_Row = 0; current_Row < cal_Bank.rows(); current_Row++) {
            if (cal_Bank.getInt("pindex", current_Row)==current_Part 
                && (cal_Bank.getInt("layer", current_Row)==1 || cal_Bank.getInt("layer", current_Row)==4 ) )  {
                pcal_plus_ecal_inner+=cal_Bank.getFloat("energy", current_Row);
            }
        }
        return 0.2 < pcal_plus_ecal_inner/p;
    }
    
    public boolean electron_z_vertex_cut(float vz) {
        // Note that target was moved between spring18 and fall18 if spring18 data is used in the future!!!
//        return vz>-13 && vz<12;
//        return vz>-6 && vz<0;
//        return vz>-13 && vz<12;
        return Math.abs(vz+3)<5;
    }
    
    public boolean pion_z_vertex_cut(float vz, double trigger_electron_vz) {
        return Math.abs(trigger_electron_vz - vz) < 20;
    }
    
    public boolean proton_z_vertex_cut(float vz, double pion_vz) {
        return Math.abs(pion_vz - vz) < 2;
    }
    
    public boolean pion_chi2pid_cut(int current_Part, HipoDataBank rec_Bank) {
        float chi2pid = rec_Bank.getFloat("chi2pid", current_Part);
        float px = rec_Bank.getFloat("px", current_Part);
        float py = rec_Bank.getFloat("py", current_Part);
        float pz = rec_Bank.getFloat("pz", current_Part);
        double p = Math.sqrt(px*px+py*py+pz*pz);
        
        int pid = rec_Bank.getInt("pid", current_Part); // slightly different cuts for pi+ and pi-
        double C;
        if (pid==211) { 
            C = 0.88;  
        } else {
            C = 0.93;
        }
        
        if (p<2.44) {
            return chi2pid < 3*C && chi2pid > -3*C;
        } else {
            return chi2pid < C * (0.00869 + 14.98587*Math.exp(-p/1.18236) + 1.81751*Math.exp(-p/4.86394))
                && chi2pid > -3*C;
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    public boolean particle_test(int current_Part, HipoDataBank rec_Bank) {
        return true
            && theta_cut(current_Part, rec_Bank) // all particles required to be detected in FD for first publications
            ;
    }
    
    public boolean electron_test(int current_Part, double p, float vz, HipoDataBank rec_Bank, HipoDataBank cal_Bank, 
            HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank, HipoDataBank cc_Bank) {
        return true
            && p > 2.0 // higher cut ultimately enforced when we cut on y < 0.8 or y < 0.75
                // this is just to speed up processing
            && forward_detector_cut(current_Part, rec_Bank)
            && calorimeter_energy_cut(current_Part, cal_Bank) 
            && calorimeter_sampling_fraction_cut(current_Part, p, cal_Bank)
            && calorimeter_diagonal_cut(current_Part, p, cal_Bank)
//            && electron_z_vertex_cut(vz)
            && pcal_fiducial_cut(current_Part, cal_Bank)
            && dc_fiducial_cut(current_Part, rec_Bank, track_Bank, traj_Bank, run_Bank)
//            && nphe_cut(current_Part, cc_Bank) // legacy cut used in the analysis note to check the effect
                ;
    }
    
    public boolean pion_test(int current_Part, int pid, float vz, double trigger_electron_vz, HipoDataBank rec_Bank, 
            HipoDataBank cal_Bank, HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        
        float px = rec_Bank.getFloat("px", current_Part);
        float py = rec_Bank.getFloat("py", current_Part);
        float pz = rec_Bank.getFloat("pz", current_Part);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
            && p > 1.00
            && p < 5 // this wasn't used in the dihadron publication but was used in the submitted single pion
            && forward_detector_cut(current_Part, rec_Bank)
//            && pion_z_vertex_cut(vz, trigger_electron_vz)
            && pion_chi2pid_cut(current_Part, rec_Bank)
            && dc_fiducial_cut(current_Part, rec_Bank, track_Bank, traj_Bank, run_Bank)
              ;
    }
    
    public boolean kaon_test(int current_Part, int pid, float vz, double trigger_electron_vz, HipoDataBank rec_Bank, 
            HipoDataBank cal_Bank, HipoDataBank track_Bank, HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        
        float px = rec_Bank.getFloat("px", current_Part);
        float py = rec_Bank.getFloat("py", current_Part);
        float pz = rec_Bank.getFloat("pz", current_Part);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
            && p > 1.00
            && p < 5 // this wasn't used in the dihadron publication but was used in the submitted single pion
            && forward_detector_cut(current_Part, rec_Bank)
//            && pion_z_vertex_cut(vz, trigger_electron_vz)
            && pion_chi2pid_cut(current_Part, rec_Bank)
            && dc_fiducial_cut(current_Part, rec_Bank, track_Bank, traj_Bank, run_Bank)
              ;
    }
    
    public boolean proton_test(int current_Part, int pid, float vz, double pion_vz, 
            HipoDataBank rec_Bank, HipoDataBank cal_Bank, HipoDataBank track_Bank, 
            HipoDataBank traj_Bank, HipoDataBank run_Bank) {
        
        float px = rec_Bank.getFloat("px", current_Part);
        float py = rec_Bank.getFloat("py", current_Part);
        float pz = rec_Bank.getFloat("pz", current_Part);
        double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
        
        return true
            && p > 0.4
//            && proton_z_vertex_cut(vz, pion_vz)
            && forward_detector_cut(current_Part, rec_Bank)
            && dc_fiducial_cut(current_Part, rec_Bank, track_Bank, traj_Bank, run_Bank)
              ;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    @Override
    public PhysicsEvent getPhysicsEvent(DataEvent event) {
        
        if (banks_test(event)) {
            PhysicsEvent physEvent = new PhysicsEvent();
            // load the hipo banks
            // assumption is we are using trains which would require all of these banks to exist
            HipoDataBank rec_Bank = (HipoDataBank) event.getBank("REC::Particle"); // load particle bank
            HipoDataBank cal_Bank = (HipoDataBank) event.getBank("REC::Calorimeter");
            HipoDataBank cc_Bank = (HipoDataBank) event.getBank("REC::Cherenkov");
            HipoDataBank track_Bank = (HipoDataBank) event.getBank("REC::Track");
            HipoDataBank traj_Bank = (HipoDataBank) event.getBank("REC::Traj");
            HipoDataBank run_Bank = (HipoDataBank) event.getBank("RUN::config");
            
            double trigger_electron_vz = -99;
            double pion_vz = -99;
            double p_max = 0;
            int p_max_index = -99; // find the index of the highest energy electorn before FD cut
            
            // cycle over the particles in recBank and investigate electron and pion IDs
            for (int current_Part = 0; current_Part < rec_Bank.rows(); current_Part++) {
                int pid = rec_Bank.getInt("pid", current_Part);
                if (pid!=11) { continue; } // find electron candidates assigned by EventBUilder
                float px = rec_Bank.getFloat("px", current_Part);
                float py = rec_Bank.getFloat("py", current_Part);
                float pz = rec_Bank.getFloat("pz", current_Part);
                double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
                float vx = rec_Bank.getFloat("vx",current_Part);
                float vy = rec_Bank.getFloat("vy",current_Part);
                float vz = rec_Bank.getFloat("vz",current_Part);
                if (p>p_max) { // searching for the highest momentum electron to use as the "trigger electron" 
                    p_max = p;
                    p_max_index = current_Part;
                    trigger_electron_vz = vz;
                }
            }
            if (p_max_index >= 0 && highest_e_in_fd_cut(p_max_index, rec_Bank)) { 
                // require that the highest momentum electron be in the forward detector 
                // THIS MAY BE MODIFIED IN A FUTURE ANALYSIS
                float px = rec_Bank.getFloat("px", p_max_index);
                float py = rec_Bank.getFloat("py", p_max_index);
                float pz = rec_Bank.getFloat("pz", p_max_index);
                double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
                float vx = rec_Bank.getFloat("vx",p_max_index);
                float vy = rec_Bank.getFloat("vy",p_max_index);
                float vz = rec_Bank.getFloat("vz",p_max_index);
                if (particle_test(p_max_index, rec_Bank)
                        && electron_test(p_max_index, p, vz, rec_Bank, cal_Bank, track_Bank, traj_Bank, 
                                run_Bank, cc_Bank)) {
                    // this checks all of the PID requirements, if it passes all of them the electron is 
                    // added to the event below
                    Particle part = new Particle(11,px,py,pz,vx,vy,vz);
                    physEvent.addParticle(part);

                }
            }

            for (int current_Part = 0; current_Part < rec_Bank.rows(); current_Part++) {
                int pid = rec_Bank.getInt("pid", current_Part);
                
                if ((Math.abs(pid)!=211) || trigger_electron_vz == -99) { continue; }
                // requires the particle to be pion by EventBuilder and for an electron to have been assigned to event
                // if no electron was assigned we just skip
                
                // load momenta and vertices
                float px = rec_Bank.getFloat("px", current_Part);
                float py = rec_Bank.getFloat("py", current_Part);
                float pz = rec_Bank.getFloat("pz", current_Part);
                double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
                float vx = rec_Bank.getFloat("vx",current_Part);
                float vy = rec_Bank.getFloat("vy",current_Part);
                float vz = rec_Bank.getFloat("vz",current_Part);
                pion_vz = vz;
                if (particle_test(current_Part, rec_Bank) 
                    && pion_test(current_Part, pid, vz, trigger_electron_vz, rec_Bank, cal_Bank, 
                    track_Bank, traj_Bank, run_Bank)) {
                    // check for pion PID
                   
                   Particle part = new Particle(pid,px,py,pz,vx,vy,vz);
                   physEvent.addParticle(part);   
               }
            }
            
            for (int current_Part = 0; current_Part < rec_Bank.rows(); current_Part++) {
                int pid = rec_Bank.getInt("pid", current_Part);
                
                if ((Math.abs(pid)!=321) || trigger_electron_vz == -99) { continue; }
                // requires the particle to be pion by EventBuilder and for an electron to have been assigned to event
                // if no electron was assigned we just skip
                
                // load momenta and vertices
                float px = rec_Bank.getFloat("px", current_Part);
                float py = rec_Bank.getFloat("py", current_Part);
                float pz = rec_Bank.getFloat("pz", current_Part);
                double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
                float vx = rec_Bank.getFloat("vx",current_Part);
                float vy = rec_Bank.getFloat("vy",current_Part);
                float vz = rec_Bank.getFloat("vz",current_Part);
                pion_vz = vz;
                if (particle_test(current_Part, rec_Bank) 
                    && kaon_test(current_Part, pid, vz, trigger_electron_vz, rec_Bank, cal_Bank, 
                    track_Bank, traj_Bank, run_Bank)) {
                    // check for pion PID
                   
                   Particle part = new Particle(pid,px,py,pz,vx,vy,vz);
                   physEvent.addParticle(part);   
               }
            }
            
            
            for (int current_Part = 0; current_Part < rec_Bank.rows(); current_Part++) {
                int pid = rec_Bank.getInt("pid", current_Part);
                
                if (pid!=2212 || pion_vz == -99) { continue; }
                
                // requires the particle to be proton by EventBuilder & for an electron to have been assigned to event
                
                // load momenta and vertices
                float px = rec_Bank.getFloat("px", current_Part);
                float py = rec_Bank.getFloat("py", current_Part);
                float pz = rec_Bank.getFloat("pz", current_Part);
                double p = Math.sqrt(Math.pow(px,2)+Math.pow(py,2)+Math.pow(pz,2));
                float vx = rec_Bank.getFloat("vx",current_Part);
                float vy = rec_Bank.getFloat("vy",current_Part);
                float vz = rec_Bank.getFloat("vz",current_Part);
                if (particle_test(current_Part, rec_Bank) 
                    && proton_test(current_Part, pid, vz, trigger_electron_vz, rec_Bank, cal_Bank, 
                    track_Bank, traj_Bank, run_Bank)) {
                   
                   Particle part = new Particle(pid,px,py,pz,vx,vy,vz);
                   physEvent.addParticle(part);   
               }
            }
            
            return physEvent;
        }
    return new PhysicsEvent(this.mybeam);
    }
}