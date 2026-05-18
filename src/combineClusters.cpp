#include "mex.hpp"
#include "mexAdapter.hpp"

#include <algorithm>
#include <vector>
#include <memory>

class MexFunction : public matlab::mex::Function {
private:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    matlab::data::ArrayFactory factory;
    
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        
        // Check for proper number of arguments
        if (inputs.size() != 3) {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({
                    factory.createScalar("three inputs required: labelmat, connmat, total")
                }));
            return;
        }
        if (outputs.size() != 1) {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({
                    factory.createScalar("one output required")
                }));
            return;
        }
        
        // Validate input types
        if (inputs[0].getType() != matlab::data::ArrayType::UINT32) {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({
                    factory.createScalar("first input must be a matrix of uint32")
                }));
            return;
        }
        
        if (inputs[1].getType() != matlab::data::ArrayType::LOGICAL) {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({
                    factory.createScalar("second input must be logical matrix")
                }));
            return;
        }
        
        if (inputs[2].getType() != matlab::data::ArrayType::UINT32 || inputs[2].getNumberOfElements() != 1) {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({
                    factory.createScalar("third input must be a scalar of uint32")
                }));
            return;
        }
        
        // Get dimensions of labelmat
        matlab::data::ArrayDimensions dims = inputs[0].getDimensions();
        size_t spatdimlength = dims[0];  // rows (channels)
        size_t timefreqlength = dims[1]; // columns (time/freq points)
        
        // Validate connmat matrix dimensions
        matlab::data::ArrayDimensions neighbourDims = inputs[1].getDimensions();
        if (neighbourDims.size() != 2 || neighbourDims[0] != spatdimlength || neighbourDims[1] != spatdimlength) {
            matlabPtr->feval(u"error", 0, 
                std::vector<matlab::data::Array>({
                    factory.createScalar("second input must be square matrix with one row and column for each channel")
                }));
            return;
        }
        
        // Get typed access to inputs (move semantics for efficiency)
        matlab::data::TypedArray<uint32_t> labelmat = std::move(inputs[0]);
        matlab::data::TypedArray<bool> connmat = std::move(inputs[1]);
        matlab::data::TypedArray<uint32_t> totalInput = std::move(inputs[2]);
        
        // Get the scalar total value (increase by 1 because indices are 1-based)
        uint32_t total = totalInput[0][0] + 1;
        
        // Create output matrix
        matlab::data::TypedArray<uint32_t> out = factory.createArray<uint32_t>(dims);
        
        // Call the computational routine
        // Create replaceby array using vector (RAII - no manual memory management)
        std::vector<uint32_t> replaceby(total);
        for (uint32_t n = 0; n < total; n++) {
            replaceby[n] = n;
        }
        
        // Iterate over channels
        for (size_t i = 0; i < spatdimlength; i++) {
            // Iterate over possible connmat for this channel
            for (size_t j = 0; j < spatdimlength; j++) {
                // Check if channels i and j are connmat
                if (connmat[i][j]) {
                    // Channel is a neighbour
                    for (size_t k = 0; k < timefreqlength; k++) {
                        uint32_t a = labelmat[i][k];
                        uint32_t b = labelmat[j][k];
                        
                        if (a > 0 && b > 0) {
                            if (replaceby[a] == replaceby[b]) {
                                continue;
                            } else if (replaceby[a] < replaceby[b]) {
                                uint32_t target = replaceby[b];
                                uint32_t replacement = replaceby[a];
                                for (uint32_t n = 0; n < total; n++) {
                                    if (replaceby[n] == target) {
                                        replaceby[n] = replacement;
                                    }
                                }
                            } else if (replaceby[b] < replaceby[a]) {
                                uint32_t target = replaceby[a];
                                uint32_t replacement = replaceby[b];
                                for (uint32_t n = 0; n < total; n++) {
                                    if (replaceby[n] == target) {
                                        replaceby[n] = replacement;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Copy and sort replaceby, retain only unique elements
        std::vector<uint32_t> replacebySorted(replaceby.begin(), replaceby.end());
        std::sort(replacebySorted.begin(), replacebySorted.end());
        auto it = std::unique(replacebySorted.begin(), replacebySorted.end());
        replacebySorted.erase(it, replacebySorted.end());
        
        // Generate sequential cluster numbers
        size_t sortedSize = replacebySorted.size();
        std::vector<uint32_t> clusternums(sortedSize);
        for (size_t n = 0; n < sortedSize; n++) {
            clusternums[n] = static_cast<uint32_t>(n); // first element will be 0 (no cluster)
        }
        
        // Fill output matrix
        for (size_t i = 0; i < spatdimlength; i++) {
            for (size_t j = 0; j < timefreqlength; j++) {
                uint32_t val = labelmat[i][j];
                if (val > 0) {
                    // Look for the cluster number
                    for (size_t n = 0; n < sortedSize; n++) {
                        if (replacebySorted[n] == replaceby[val]) {
                            out[i][j] = clusternums[n];
                            break;
                        }
                    }
                }
            }
        }

        // Return output
        outputs[0] = std::move(out);
    }
};
