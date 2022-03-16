#include "FileUtil.h"

std::string FileUtil::getFileString(std::string path,bool withEnd) {
    FILE* fp = 0;
    fopen_s(&fp, path.c_str(), "rb");
    if (fp != nullptr)
    {
        fseek(fp, 0, SEEK_END);
        int fileSize = ftell(fp);
        fseek(fp, 0, SEEK_SET);
        char* buffer;
        if (withEnd) {
            buffer= (char*)malloc(fileSize + 1);
            buffer[fileSize] = '\0';
        }
        else {
            buffer = (char*)malloc(fileSize);
        }     
        
        int readsize = fread(buffer, 1, fileSize, fp);
        if (readsize < fileSize)
        {
            printf("ÎÄ¼þ¶ÁÈ¡´íÎó\n");
            return "";
        }
        fclose(fp);
        std::string str(buffer);
        free(buffer);
        return str;

    }
    return "";
}

std::string FileUtil::getRootPath() {
    char buff[261];
    DWORD nsize = sizeof(buff) / sizeof(char);
    DWORD n = GetModuleFileNameA(NULL, buff, nsize);
    std::string exePath(buff);
    int pos = exePath.find_first_of("x64");
    std::string root = exePath.substr(0, pos);
    return root;
}