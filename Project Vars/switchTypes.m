% Helper function to isolate the switch statement for constructing
% variables.
function [kelv ross force refl rere] = switchTypes(option)

    switch option
        case 1  %All forcing
                kelv = true;
                ross = true;

                force = true;
                refl = true;
                rere = true;
                display('All Forcing');
        case 2  %Kelvin Forced
                kelv = true;
                ross = false;

                force = true;
                refl = false;
                rere = false;
                display('Kelvin - Forced');
        case 3  %Kelvin Reflected
                kelv = true;
                ross = false;

                force = false;
                refl = true;
                rere = false;
                display('Kelvin - Reflected');

        case 4  %Kelvin Reflected twice
                kelv = true;
                ross = false;

                force = false;
                refl = false;
                rere = true;                
                display('Kelvin - Reflected x 2');

        case 5  %Rossby Forced
                kelv = false;
                ross = true;

                force = true;
                refl = false;
                rere = false;                
                display('Rossby - Forced');

        case 6  %Rossby Reflected
                kelv = false;
                ross = true;

                force = false;
                refl = true;
                rere = false;                
                display('Rossby - Reflected');

        case 7  %Rossby Reflected twice
                kelv = false;
                ross = true;

                force = false;
                refl = false;
                rere = true;                
                display('Rossby - Reflected x 2');

    end
end

