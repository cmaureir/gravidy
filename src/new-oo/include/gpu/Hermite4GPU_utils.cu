#include "Hermite4GPU.cuh"

void Hermite4GPU::get_kernel_error(){
    #ifdef KERNEL_ERROR_DEBUG
        std::cerr << "[Error] : ";
        std::cerr << cudaGetErrorString(cudaGetLastError()) << std::endl;
    #endif
}

void Hermite4GPU::gpu_timer_start(){
    cudaEventRecord(start);
}

float HermiteGPU::gpu_timer_stop(string f){
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float msec = 0;
    cudaEventElapsedTime(&msec, start, stop);
    #if KERNEL_TIME
    if (f != "")
        std::cout << "[Time] " << f << " : " << msec << " msec" << std::endl;
    #endif
    return msec;
}
