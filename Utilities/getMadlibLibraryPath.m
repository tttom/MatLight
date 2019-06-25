function madlibLibraryPath=getMadlibLibraryPath()
    basePath=getLibraryPath();
    
    switch (upper(computer()))
        case 'PCWIN64'
            madlibLibraryPath=strcat(basePath,'Madlib\');
        otherwise
            %Assume 32 bit Windows
            madlibLibraryPath=strcat(basePath,'Madlib32\');
    end
end