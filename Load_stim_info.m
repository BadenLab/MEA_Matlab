function [Epochs] = Load_stim_info (Epochs, stimulus_path)

epochs = length(Epochs.Epoch_code);
cd(stimulus_path);

compare = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7];
compare = single(compare);

for ii = 1:epochs
    
    switch Epochs.Epoch_code(ii)
        
    case compare(1)
        
        load('FFF430');
        Epochs.FFF430 = FFF430;
        Epochs.Epochtext(1,ii) = {'FFF430'};
        
    case compare(2)
       load('FFF480');
       Epochs.FFF480 = FFF480;
       Epochs.Epochtext(1,ii) = {'FFF480'};
       
    case compare(3)
        load('FFF505');
        Epochs.FFF505 = FFF505;
        Epochs.Epochtext(1,ii) = {'FFF505'};
        
    case compare(4)
      load('FFF560');
      Epochs.FFF560 = FFF560;
      Epochs.Epochtext(1,ii) = {'FFF560'};
      
    case compare(5)
        load('FFF630');
        Epochs.FFF630 = FFF630;
        Epochs.Epochtext(1,ii) = {'FFF630'};
        
    case compare(6)
        load('FFF');
        Epochs.FFF = FFF;
        Epochs.Epochtext(1,ii) = {'FFF'};
        
    case compare(7)
        load('Noise380');
        Epochs.Noise380 = Noise380;
        load('Noise430');
        Epochs.Noise430 = Noise430;
        load('Noise480');
        Epochs.Noise480 = Noise480;
        load('Noise505');
        Epochs.Noise505 = Noise505;
        load('Noise560');
        Epochs.Noise560 = Noise560;
        load('Noise630');
        Epochs.Noise630 = Noise630;
        Epochs.Epochtext(1,ii) = {'CNoise'};
        
        
    case compare(8)
        load('chirp1');
        Epochs.Chirp = chirp1;
        Epochs.Epochtext(1,ii) = {'Chirp'};
        
    case compare(9)
        load('FFF380');
        Epochs.FFF380 = FFF380;
        Epochs.Epochtext(1,ii) = {'FFF380'};
        
    case compare(10)
        load('SS430');
        Epochs.SS430 = SS430;
        Epochs.Epochtext(1,ii) = {'SS430'};
        
    case compare(11)
        load('SS480');
        Epochs.SS480 = SS480;
        Epochs.Epochtext(1,ii) = {'SS480'};
        
    case compare(12)
        load('SS505');
        Epochs.SS505 = SS505;
        Epochs.Epochtext(1,ii) = {'SS505'};
        
    case compare(13)
        load('SS560');
        Epochs.SS560 = SS560;
        Epochs.Epochtext(1,ii) = {'SS560'};
        
    case compare(14)
        load('SS630');
        Epochs.SS630 = SS630;
        Epochs.Epochtext(1,ii) = {'SS630'};
        
            
            
    end
end







end