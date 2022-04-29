//-------------------------------------------------------------------------------------
// DirectXMatrixStack.h -- DirectXMath C++ Matrix Stack
//
// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.
//
// http://go.microsoft.com/fwlink/?LinkID=615560
//-------------------------------------------------------------------------------------

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <new>

#ifdef _WIN32
#include <malloc.h>
#endif

#include <DirectXMath.h>


namespace DirectX
{
    class MatrixStack
    {
    public:
        MatrixStack(size_t startSize = 16) noexcept(false) :
            m_stackSize(0),
            m_current(0),
            m_stack(nullptr)
        {
            assert(startSize > 0);
            Allocate(startSize);
            LoadIdentity();
        }

        MatrixStack(MatrixStack&&) = default;
        MatrixStack& operator= (MatrixStack&&) = default;

        MatrixStack(MatrixStack const&) = delete;
        MatrixStack& operator= (MatrixStack const&) = delete;

        const XMMATRIX XM_CALLCONV Top() const noexcept { return m_stack[m_current]; }
        const XMMATRIX* GetTop() const noexcept { return &m_stack[m_current]; }

        size_t Size() const noexcept { return (m_current + 1); }

        void Pop()
        {
            if (m_current > 0)
            {
                --m_current;
            }
        }

        void Push()
        {
            ++m_current;

            if (m_current >= m_stackSize)
            {
                Allocate(m_stackSize * 2);
            }

            // Replicate the original top of the matrix stack.
            m_stack[m_current] = m_stack[m_current - 1];
        }

        // Loads identity into the top of the matrix stack.
        void LoadIdentity() noexcept
        {
            m_stack[m_current] = XMMatrixIdentity();
        }

        // Load a matrix into the top of the matrix stack.
        void XM_CALLCONV LoadMatrix(FXMMATRIX matrix) noexcept
        {
            m_stack[m_current] = matrix;
        }

        // Multiply a matrix by the top of the stack, store result in top.
        void XM_CALLCONV MultiplyMatrix(FXMMATRIX matrix) noexcept
        {
            m_stack[m_current] = XMMatrixMultiply(m_stack[m_current], matrix);
        }

        // Pre-multiplies a matrix by the top of the stack, store result in top.
        void XM_CALLCONV MultiplyMatrixLocal(FXMMATRIX matrix) noexcept
        {
            m_stack[m_current] = XMMatrixMultiply(matrix, m_stack[m_current]);
        }

        // Add a rotation about X to stack top.
        void XM_CALLCONV RotateX(float angle) noexcept
        {
            XMMATRIX mat = XMMatrixRotationX(angle);
            m_stack[m_current] = XMMatrixMultiply(m_stack[m_current], mat);
        }

        void XM_CALLCONV RotateXLocal(float angle) noexcept
        {
            XMMATRIX mat = XMMatrixRotationX(angle);
            m_stack[m_current] = XMMatrixMultiply(mat, m_stack[m_current]);
        }

        // Add a rotation about Y to stack top.
        void XM_CALLCONV RotateY(float angle) noexcept
        {
            XMMATRIX mat = XMMatrixRotationY(angle);
            m_stack[m_current] = XMMatrixMultiply(m_stack[m_current], mat);
        }

        void XM_CALLCONV RotateYLocal(float angle) noexcept
        {
            XMMATRIX mat = XMMatrixRotationY(angle);
            m_stack[m_current] = XMMatrixMultiply(mat, m_stack[m_current]);
        }

        // Add a rotation about Z to stack top.
        void XM_CALLCONV RotateZ(float angle) noexcept
        {
            XMMATRIX mat = XMMatrixRotationZ(angle);
            m_stack[m_current] = XMMatrixMultiply(m_stack[m_current], mat);
        }

        void XM_CALLCONV RotateZLocal(float angle) noexcept
        {
            XMMATRIX mat = XMMatrixRotationZ(angle);
            m_stack[m_current] = XMMatrixMultiply(mat, m_stack[m_current]);
        }

        // Add a rotation around an axis to stack top.
        void XM_CALLCONV RotateAxis(FXMVECTOR axis, float angle) noexcept
        {
            XMMATRIX mat = XMMatrixRotationAxis(axis, angle);
            m_stack[m_current] = XMMatrixMultiply(m_stack[m_current], mat);
        }

        void XM_CALLCONV RotateAxisLocal(FXMVECTOR axis, float angle) noexcept
        {
            XMMATRIX mat = XMMatrixRotationAxis(axis, angle);
            m_stack[m_current] = XMMatrixMultiply(mat, m_stack[m_current]);
        }

        // Add a rotation by roll/pitch/yaw to the stack top.
        void RotateRollPitchYaw(float pitch, float yaw, float roll) noexcept
        {
            XMMATRIX mat = XMMatrixRotationRollPitchYaw(pitch, yaw, roll);
            m_stack[m_current] = XMMatrixMultiply(m_stack[m_current], mat);
        }

        void RotateRollPitchYawLocal(float pitch, float yaw, float roll) noexcept
        {
            XMMATRIX mat = XMMatrixRotationRollPitchYaw(pitch, yaw, roll);
            m_stack[m_current] = XMMatrixMultiply(mat, m_stack[m_current]);
        }

        // Add a rotation by a quaternion stack top.
        void XM_CALLCONV RotateByQuaternion(FXMVECTOR quat) noexcept
        {
            XMMATRIX mat = XMMatrixRotationQuaternion(quat);
            m_stack[m_current] = XMMatrixMultiply(m_stack[m_current], mat);
        }

        void XM_CALLCONV RotateByQuaternionLocal(FXMVECTOR quat) noexcept
        {
            XMMATRIX mat = XMMatrixRotationQuaternion(quat);
            m_stack[m_current] = XMMatrixMultiply(mat, m_stack[m_current]);
        }

        // Add a scale to the stack top.
        void Scale(float x, float y, float z) noexcept
        {
            XMMATRIX mat = XMMatrixScaling(x, y, z);
            m_stack[m_current] = XMMatrixMultiply(m_stack[m_current], mat);
        }

        void ScaleLocal(float x, float y, float z) noexcept
        {
            XMMATRIX mat = XMMatrixScaling(x, y, z);
            m_stack[m_current] = XMMatrixMultiply(mat, m_stack[m_current]);
        }

        // Add a translation to the stack top.
        void Translate(float x, float y, float z) noexcept
        {
            XMMATRIX mat = XMMatrixTranslation(x, y, z);
            m_stack[m_current] = XMMatrixMultiply(m_stack[m_current], mat);
        }

        void TranslateLocal(float x, float y, float z) noexcept
        {
            XMMATRIX mat = XMMatrixTranslation(x, y, z);
            m_stack[m_current] = XMMatrixMultiply(mat, m_stack[m_current]);
        }

    private:

        struct matrix_deleter
        {
            void operator()(void* p) noexcept
            {
#ifdef _WIN32
                _aligned_free(p);
#else
                free(p);
#endif
            }
        };

        void Allocate(size_t newSize)
        {
#ifdef _WIN32
            void* ptr = _aligned_malloc(newSize * sizeof(XMMATRIX), 16);
#else
            // This C++17 Standard Library function is currently NOT
            // implemented for the Microsoft Standard C++ Library.
            void* ptr = aligned_alloc(16, newSize * sizeof(XMMATRIX));
#endif
            if (!ptr)
                throw std::bad_alloc();

            if (m_stack)
            {
                assert(newSize >= m_stackSize);
                memcpy(ptr, m_stack.get(), sizeof(XMMATRIX) * m_stackSize);
            }

            m_stack.reset(reinterpret_cast<XMMATRIX*>(ptr));
            m_stackSize = newSize;
        }

        size_t										m_stackSize;
        size_t										m_current;
        std::unique_ptr<XMMATRIX[], matrix_deleter>	m_stack;
    };
} // namespace DirectX
